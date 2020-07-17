#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec hermite(const arma::vec& x, int n) {
  if (n == 0) {
    return arma::ones(x.n_elem);
  } else if (n == 1) {
    return x;
  } else {
    return x % hermite(x, n - 1) - (n - 1) * hermite(x, n - 2);
  }
}
