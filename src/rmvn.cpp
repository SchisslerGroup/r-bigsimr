#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat rmvn(int n, arma::colvec mu, arma::mat sigma) {
  return arma::mvnrnd(mu, sigma, n).t();
}
