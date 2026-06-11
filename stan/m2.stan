// adapted from:
// https://gptools-stan.readthedocs.io/docs/logistic_regression/logistic_regression.html

// from m1: add padding points

functions {
    #include gptools/util.stan
    #include gptools/fft1.stan
}

data {
    // int n, p;
    // NO padding. We WANT the periodic boundary conditions here.
    int n;
    int n_fill;
    array [n_fill] int<lower=0> n_i; // trials at week i
    array [n] int<lower=0> obs_idx; // indices in n_i with real observations
    array [n_fill] int<lower=0> y;
    real<lower=0> epsilon;
    // TODO: include year-level effect
}

parameters {
    real<lower=0> sigma;
    real<lower=log(2), upper=log(n_fill / 2.0)> log_length_scale;
    vector[n_fill] raw;
}

transformed parameters {
    real length_scale = exp(log_length_scale);
    // Evaluate the covariance kernel in the Fourier domain. We add the nugget
    // variance because the Fourier transform of a delta function is a constant.
    vector[n_fill %/% 2 + 1] cov_rfft = gp_periodic_exp_quad_cov_rfft(
        n_fill, sigma, length_scale, n_fill) + epsilon;
    // Transform the "raw" parameters to the latent log odds ratio.
    vector[n_fill] z = gp_inv_rfft(raw, zeros_vector(n_fill), cov_rfft);
}

model {
    raw ~ std_normal();
    y[obs_idx] ~ binomial_logit(n_i[obs_idx], z[obs_idx]);
    sigma ~ normal(0, 3);
}

generated quantities {
    vector[n_fill] proba = inv_logit(z);
}
