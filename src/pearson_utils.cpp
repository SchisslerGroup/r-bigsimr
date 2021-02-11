#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double Hn(double x, int n) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  }

  double Hkp1 = 0.0;
  double Hk = x;
  double Hkm1 = 1.0;
  for (int k = 2; k <= n; k++) {
    Hkp1 = x*Hk - (k-1) * Hkm1;
    Hkm1 = Hk;
    Hk = Hkp1;
  }

  return Hkp1;
}

// [[Rcpp::export]]
double Hp(double x, int n) {
  if (std::isinf(x)) {
    return 0.0;
  } else {
    NumericVector xR {x};
    double dnorm_x = dnorm(xR)[0];
    return Hn(x, n) * dnorm_x;
  }
}

// [[Rcpp::export]]
double Gn0d(const arma::vec& A, const arma::vec& B,
            const arma::vec& a, const arma::vec& b,
            double sAsB_inv, int k) {

  if (k == 0) {
    return 0.0;
  }

  int M = A.n_elem;
  int N = B.n_elem;

  double accu = 0.0;
  double r11, r00, r01, r10;

  for (int r = 0; r < M; r++) {
    for (int s = 0; s < N; s++) {
      r11 = Hp(a[r+1], k-1) * Hp(b[s+1], k-1);
      r00 = Hp(a[r],   k-1) * Hp(b[s],   k-1);
      r01 = Hp(a[r],   k-1) * Hp(b[s+1], k-1);
      r10 = Hp(a[r+1], k-1) * Hp(b[s],   k-1);
      accu += A[r] * B[s] * (r11 + r00 - r01 - r10);
    }
  }
  return accu * sAsB_inv;
}
