#include <RcppGSL.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

// see pg. 141+ of the GSL docs pdf for solving linear systems
// cholesky on p 154
// gsl_sf_lnchoose_e on pg 60

// declare a dependency on the RcppGSL package; also activates plugin
// (but not needed when ’LinkingTo: RcppGSL’ is used with a package)
//
// [[Rcpp::depends(RcppGSL)]]
// tell Rcpp to turn this into a callable function called ’fastLm’
//
// [[Rcpp::export]]
Rcpp::NumericVector lpost_opt(const RcppGSL::Vector & pvec, const RcppGSL::Vector & y) {
  // int n = X.nrow(), k = X.ncol(); // row and column dimension
  
  //pull out parameters
  
  //get L
  
  // log norm week effects
  
  // log norm yrweek effects
  
  // add those two to other priors on l, wk_sigma, and intercept
  
  // get linear predictor: wkv[wk_i] + yrwkv + intercept
  
  // get binomial coefficients
  
  // get likelihood: sum(lc + k*log(p) + (n-k) * log(q))
  
  // double chisq; // assigned but not returned
  // RcppGSL::Vector coef(k); // to hold the coefficient vector
  // RcppGSL::Matrix cov(k,k); // and the covariance matrix
  // // the actual fit requires working memory which we allocate and then free
  // gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, k);
  // gsl_multifit_linear (X, y, coef, cov, &chisq, work);
  // gsl_multifit_linear_free (work);
  // // assign diagonal to a vector, then take square roots to get std.error
  // Rcpp::NumericVector std_err;
  // std_err = gsl_matrix_diagonal(cov); // need two step decl. and assignment
  // std_err = Rcpp::sqrt(std_err); // sqrt() is an Rcpp sugar function
  Rcpp::NumericVector res(5);
  return res;
}