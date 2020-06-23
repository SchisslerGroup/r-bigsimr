#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat rmvuu(int n, arma::mat rho) {
  arma::mat Z = arma::mvnrnd(arma::zeros(rho.n_cols), rho, n).t();
  return arma::normcdf(Z);
}


// [[Rcpp::export]]
arma::mat rmvu(int n, arma::mat rho, arma::rowvec min, arma::rowvec max) {
  arma::mat l = arma::ones(n, 1) * min;
  arma::mat u = arma::ones(n, 1) * max;

  return rmvuu(n, rho) % (u - l) + l;
}
