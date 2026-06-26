functions {
  matrix getL(real ell, real wk_sigma) {
    matrix[52,52] cov_mat = rep_matrix(0, 52, 52);
    // cholesky_factor_cov[52] L;
    matrix[52,52] L = rep_matrix(0, 52, 52);
    
    for (i in 2:27) {
      cov_mat[i, 1] = -(i-1)^2 / (2*ell^2);
    }
    
    cov_mat[27:52, 1] = cov_mat[reverse(linspaced_int_array(26, 2, 27)), 1];
    
    vector[52] cv = cov_mat[,1];
    
    array[52] int rv = linspaced_int_array(52, 0, 51);
    rv[1] = 52;
    
    for (j in 2:52) {
      // cv = cv[rv];
      cov_mat[,j] = cov_mat[,j-1][rv];
    }
   
    cov_mat = wk_sigma^2 * exp(cov_mat); 
    
    // print(cov_mat[1:5,1:5]);
    L = cholesky_decompose(cov_mat);
    
    return L;
  }
}

data {
  
  int<lower=0> N;
  array[N] int<lower=0> k;
  array[N] int<lower=0> n;
  array[N] int<lower=1, upper=52> wk_i;
  
}

parameters {
  real<lower=0> ell;
  real<lower=0> wk_sigma;
  real intercept;
  
  vector[52] wkv;
  vector[N] yrwkv;
  
}

transformed parameters {
  // cholesky_factor_cov[52] L = getL(ell, wk_sigma);
  // vector[52] z = mdivide_left_tri_low(L, wkv);
  // 
  // vector[N] lin_pred = intercept + wkv[wk_i] + yrwkv;
}

model {

  
  matrix[52,52] L = getL(ell, wk_sigma);
  vector[52] z = mdivide_left_tri_low(L, wkv);
  
  vector[N] lin_pred = intercept + wkv[wk_i] + yrwkv;
  
  z ~ std_normal();
  
  yrwkv ~ normal(0, .3);
  
  ell ~ inv_gamma(2,2);
  
  wk_sigma ~ exponential(10);
  
  intercept ~ normal(0,2);
  
  k ~ binomial_logit(n, lin_pred);
  
}
