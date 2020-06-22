#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#define PEARSON_SPEARMAN 0
#define PEARSON_KENDALL  1
#define SPEARMAN_PEARSON 2
#define SPEARMAN_KENDALL 3
#define KENDALL_PEARSON  4
#define KENDALL_SPEARMAN 5

// [[Rcpp::export]]
arma::mat CXX_cor2cor(const arma::mat& X, int CASE) {

  arma::mat rho(size(X), arma::fill::zeros);

  switch (CASE) {
    case PEARSON_SPEARMAN:
      rho += (6 / M_PI) * arma::asin(X / 2);
      break;
    case PEARSON_KENDALL:
      rho += (2 / M_PI) * arma::asin(X);
      break;
    case SPEARMAN_PEARSON:
      rho += 2 * arma::sin(X * M_PI / 6);
      break;
    case SPEARMAN_KENDALL:
      rho += (2 / M_PI) * arma::asin(2 * arma::sin(X * M_PI / 6));
      break;
    case KENDALL_PEARSON:
      rho += arma::sin(X * M_PI / 2);
      break;
    case KENDALL_SPEARMAN:
      rho += (6 / M_PI) * arma::asin(arma::sin(X * M_PI / 2) / 2);
      break;
    default:
      rho += X;
  }

  return rho;
}
