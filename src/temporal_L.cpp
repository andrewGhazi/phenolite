#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
mat getL(double l , double wk_sigma) {
  
  mat cov_mat = zeros(52, 52);
  vec cv = zeros(52);
  double linv = 1/(2*pow(l, 2.0));
  
  int j = 51;
  for (int i = 1; i<28; ++i) {
    // ^ i starting at 1 here is not an indexing error, really mean to start at the second element.
    
    cv[i] = -pow(i, 2.0) * linv;
    cv[j] = cv[i];
    j -= 1;
  }
  
  for (int i = 0; i<52; ++i) {
    cov_mat.col(i) = cv;
    
    cv = circshift(cv, 1);
  }
  
  cov_mat = pow(wk_sigma, 2.0) * exp(cov_mat); 
  
  mat L = chol(cov_mat, "lower");
  // Rcpp::Rcout << cov_mat << "\n";
  
  return L;
}

// [[Rcpp::export]]
double dexp_l(double x, double rt){
  return log(rt) - rt*x;
}

// [[Rcpp::export]]
arma::vec inv_logit(arma::vec x) {
  arma::vec res = 1 / (1 + exp(-x));
  return res;
}


Rcpp::NumericVector inv_logit(Rcpp::NumericVector x) {
  Rcpp::NumericVector res = 1 / (1 + exp(-x));
  return res;
}

double inv_logit(double x) {
  double res = 1 / (1 + exp(-x));
  return res;
}

// [[Rcpp::export]]
double lpost(Rcpp::NumericVector pvec,
             int N, 
             Rcpp::NumericVector k, 
             Rcpp::NumericVector n, 
             Rcpp::IntegerVector wk_i) {
  
  double l = pvec[0];
  double wk_sigma = pvec[1];
  double intercept = pvec[2];
  
  int n_hp = 3;
  
  Rcpp::NumericVector wkv = pvec[Rcpp::Range(3, 54)]; //.subvec(3, 54);
  // Rcpp::NumericVector yrwkv = tail(pvec,N);

  mat L = getL(l, wk_sigma);
  vec Ls = solve(L, Rcpp::as<arma::vec>(Rcpp::wrap(wkv)));

  vec wk_ll = log_normpdf(Ls);
  double wk_lls = sum(wk_ll);

  double yrwk_ll = sum(Rcpp::dnorm(tail(pvec,N), 0, 0.5, true));
  // Rcpp::NumericVector yrwk_llv = sum(Rcpp::dnorm(tail(pvec,N), 0, 0.5, true));
  // double yrwk_ll = sum(yrwk_llv);
  
  double lp = wk_lls +
              dexp_l(wk_sigma, 10) +
              dexp_l(l, 0.5) +
              yrwk_ll +
              log_normpdf(intercept);
  
  
  Rcpp::NumericVector lc = log(Rcpp::choose(n, k));
  // Rcpp::NumericVector lc = log(coefs);
  Rcpp::NumericVector x = wkv[wk_i];
  
  x = x + tail(pvec,N) + intercept;
  
  Rcpp::NumericVector p = inv_logit(x);
  Rcpp::NumericVector nklq = (n-k) * log(1-p);
  Rcpp::NumericVector klp = k*log(p);
  // Rcpp::NumericVector nklq = (n-k) * log(q);
  
  double ll = sum(lc + klp + nklq);
  
  return lp + ll;
}


// [[Rcpp::export]]
double lpost2(vec pvec,
              int N, 
              uvec k, 
              uvec n, 
              uvec wk_i) {
  
  double l = pvec[0];
  double wk_sigma = pvec[1];
  double intercept = pvec[2];
  
  int n_hp = 3;
  
  vec wkv = pvec.subvec(3, 54);

  // 52 weeks follows temporal GP with circular boundary conditions.
  mat L = getL(l, wk_sigma);
  vec Lsol = solve(L, wkv);
  double wk_lls = accu(log_normpdf(Lsol));
  // Lz = v induces desired covariance structure on V from std normally
  // distributed z. So this std normal prior on Lsol ensures wkv (the 52 weekly
  // points reflecting the average behavior of that week) has the temporal
  // correlation + joined boundaries.

  // Each observation (a particular week in a given year) gets a little bit of
  // wiggle room
  double yrwk_ll = accu(log_normpdf(pvec.tail(N) , 0, 0.4));
  
  double lp = wk_lls +
              dexp_l(wk_sigma, 10) +
              dexp_l(l, 0.4) +
              yrwk_ll +
              log_normpdf(intercept);

  // Need to go to a NumVec and back because Armadillo doesn't have lchoose
  Rcpp::NumericVector lc = Rcpp::lchoose(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(n)),
                                         Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(k)));
  
  vec lca = Rcpp::as<vec>(Rcpp::wrap(lc));
  
  vec p = inv_logit(wkv(wk_i) + pvec.tail(N) + intercept);
  vec klp = k % log(p);
  vec nklq = (n-k) % log(1-p);
  
  double ll = sum(lca + klp + nklq);
   
  return lp + ll;
}
