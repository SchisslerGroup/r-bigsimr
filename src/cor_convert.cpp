#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#define PEARSON_SPEARMAN 0
#define PEARSON_KENDALL  1
#define SPEARMAN_PEARSON 2
#define SPEARMAN_KENDALL 3
#define KENDALL_PEARSON  4
#define KENDALL_SPEARMAN 5


// [[Rcpp::export(.cor_convert_double)]]
double cor_convert_double(double X, int CASE) {
  switch (CASE) {
  case PEARSON_SPEARMAN:
    return (6 / M_PI) * asin(X / 2);
  case PEARSON_KENDALL:
    return (2 / M_PI) * asin(X);
  case SPEARMAN_PEARSON:
    return 2 * sin(X * M_PI / 6);
  case SPEARMAN_KENDALL:
    return (2 / M_PI) * asin(2 * sin(X * M_PI / 6));
  case KENDALL_PEARSON:
    return sin(X * M_PI / 2);
  case KENDALL_SPEARMAN:
    return (6 / M_PI) * asin(sin(X * M_PI / 2) / 2);
  default:
    return X;
  }
}


// [[Rcpp::export(.cor_convert_matrix)]]
arma::mat cor_convert_matrix(const arma::mat& X, int CASE) {
  switch (CASE) {
  case PEARSON_SPEARMAN:
    return (6 / M_PI) * arma::asin(X / 2);
  case PEARSON_KENDALL:
    return (2 / M_PI) * arma::asin(X);
  case SPEARMAN_PEARSON:
    return 2 * arma::sin(X * M_PI / 6);
  case SPEARMAN_KENDALL:
    return (2 / M_PI) * arma::asin(2 * arma::sin(X * M_PI / 6));
  case KENDALL_PEARSON:
    return arma::sin(X * M_PI / 2);
  case KENDALL_SPEARMAN:
    return (6 / M_PI) * arma::asin(arma::sin(X * M_PI / 2) / 2);
  default:
    return X;
  }
}
