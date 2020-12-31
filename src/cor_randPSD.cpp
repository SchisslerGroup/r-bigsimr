#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export(.cor_randPSD)]]
arma::mat cor_randPSD(int d, int k) {
  if(d == 1) return arma::ones<arma::mat>(1,1);

  arma::mat W = arma::randn(d, d);
  arma::mat S = W * W.t() + arma::diagmat(arma::randu<arma::vec>(d));
  arma::mat S2 = arma::diagmat(1 / arma::sqrt(S.diag()));

  arma::mat R = arma::clamp(S2 * S * S2, -1.0, 1.0);

  for (unsigned int i = 0; i < R.n_cols; ++i) {
    R(i, i) = 1.0;
  }
  return R;
}
