// adapted from:
// https://gptools-stan.readthedocs.io/docs/logistic_regression/logistic_regression.html

functions {
    #include gptools/util.stan
    #include gptools/fft1.stan
}

data {
    // int n, p;
    // NO padding. We WANT the periodic boundary conditions here.
    int n;
    array [n] int<lower=0> n_i; // trials at week i
    array [n] int<lower=0> y;
    real<lower=0> epsilon;
    // TODO: include year-level effect
}

parameters {
    real<lower=0> sigma;
    real<lower=log(2), upper=log(n / 2.0)> log_length_scale;
    vector[n] raw;
}

transformed parameters {
    real length_scale = exp(log_length_scale);
    // Evaluate the covariance kernel in the Fourier domain. We add the nugget
    // variance because the Fourier transform of a delta function is a constant.
    vector[n %/% 2 + 1] cov_rfft = gp_periodic_exp_quad_cov_rfft(
        n, sigma, length_scale, n) + epsilon;
    // Transform the "raw" parameters to the latent log odds ratio.
    vector[n] z = gp_inv_rfft(raw, zeros_vector(n), cov_rfft);
}

model {
    raw ~ std_normal();
    y ~ binomial_logit(n_i, z);
    sigma ~ normal(0, 3);
}

generated quantities {
    vector[n] proba = inv_logit(z);
}
