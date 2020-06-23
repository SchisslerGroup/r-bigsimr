#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Code adapted from user amoeba on stats.stackexchange.com
// from https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor


// [[Rcpp::export]]
arma::mat rcor(int d, int k = 1) {
  arma::mat W = arma::randn(d, k);
  arma::mat S = W * W.t() + arma::diagmat( arma::randu(d) );
  arma::mat S_diag_inv = arma::diagmat( 1 / arma::sqrt(S.diag()) );
  return S_diag_inv * S * S_diag_inv;
}
