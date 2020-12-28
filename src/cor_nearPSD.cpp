#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#define PERTURBATION 1.0e-9


void Jacobian(
    const arma::vec& x,
    const arma::mat& Omega12,
    const arma::mat& P,
    arma::vec& Ax) {

  int n = P.n_cols;
  int r = Omega12.n_rows;
  int s = Omega12.n_cols;

  if (r == 0) {
    // TODO: this should probably throw an exception
    Ax.zeros();
  } else if (r == n) {
    Ax = x * (1.0 + PERTURBATION);
  } else if (r < s) {
    arma::mat P1 = P.head_cols(r);
    arma::mat P2 = P.tail_cols(s);

    arma::mat H1 = arma::diagmat(x) * P1;
    arma::mat Omega = Omega12 % (H1.t() * P2);

    arma::mat HT(n, n, arma::fill::zeros);
    HT.head_cols(r) += P1 * P1.t() * H1 + P2 * Omega.t();
    HT.tail_cols(s) += P1 * Omega;

    Ax = arma::sum(P % HT, 1) + x * PERTURBATION;
  } else { // (r >= s)
    arma::mat P1 = P.head_cols(r);
    arma::mat P2 = P.tail_cols(s);

    arma::mat H2 = arma::diagmat(x) * P2;
    arma::mat Omega = (1.0 - Omega12) % (P1.t() * H2);

    arma::mat HT(n, n, arma::fill::zeros);
    HT.head_cols(r) += P2 * Omega.t();
    HT.tail_cols(s) += P2 * H2.t() * P2 + P1 * Omega;

    Ax = x * (1.0 + PERTURBATION) - arma::sum(P % HT, 1);
  }
  return;
}


void gradient(
    const arma::vec& y,
    const arma::vec& lambda,
    arma::mat P,
    const arma::vec& b0,
    double& f,
    arma::vec& Fy) {

  int n = P.n_cols;
  int r = arma::accu(lambda > 0); // number of positive eigenvalues

  if (r == 0) {
    Fy.zeros();
    f = 0;
  } else {
    // work only with the positive eigenvalues (lambda)
    arma::vec l1(lambda);
    l1.elem(find(lambda < 0)).fill(0);

    // column j is multiplied by lambda[j]
    for (int j = 0; j < n; j++) {
      P.col(j) *= std::sqrt(l1(j));
    }

    Fy = arma::sum(P % P, 1);
    f = 0.5 * arma::accu(arma::pow(l1, 2.0)) - arma::accu(b0 % y);
  }
  return;
}


void omega_mat(const arma::vec& lambda, arma::mat& Omega12) {

  int n = lambda.n_elem;
  int r = arma::accu(lambda > 0);
  int s = n - r;

  if (r == 0) {
    Omega12.set_size(0,0);
  } else if (r == n) {
    Omega12.set_size(n, n);
    Omega12.ones();
  } else { // 0 < r < n
    Omega12.set_size(r, s);
    arma::vec lambda_s = lambda.tail(s);  // non-positive eigenvalues
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < s; j++) {
        Omega12(i, j) = lambda(i) / (lambda(i) - lambda_s(j));
      }
    }
  }
  return;
}


void precond_matrix(const arma::mat& Omega12, const arma::mat& P, arma::vec& c) {

  int n = P.n_cols;
  int r = Omega12.n_rows;
  int s = Omega12.n_cols;

  c.ones();

  if (r == 0 || r == n)  {
    // do not modify anything
    return;
  }

  arma::mat H = (P % P).t();
  arma::mat H1 = H.head_rows(r);
  arma::mat H2 = H.tail_rows(s);

  if (r < s) {
    arma::mat H12 = H1.t() * Omega12;
    c  = arma::pow(arma::sum(H1, 0).t(), 2.0);
    c += 2.0 * arma::sum(H12 % H2.t(), 1);
  } else {
    arma::mat H12 = (1 - Omega12) * H2;
    c  = arma::pow(arma::sum(H, 0).t(), 2.0);
    c -= arma::pow(arma::sum(H2, 0).t(), 2.0);
    c -= 2 * arma::sum((H1 % H12), 0).t();
  }
  c.elem(find(c < 1.0e-8)).fill(1.0e-8);
  return;
}


void pre_cg(const arma::vec& b, const arma::vec& c, const arma::mat& Omega12,
            const arma::mat& P, const double& precg_err_tol, const int& maxit, arma::vec& p) {

  int n = P.n_cols;

  double tolb = precg_err_tol * arma::norm(b, 2);

  arma::vec rv(b);
  arma::vec z = rv / c;
  double rz1 = arma::accu(rv % z);
  double rz2 = 1.0;

  p.zeros();
  arma::vec d(z);
  arma::vec w(n);

  for (int k = 1; k <= maxit; ++k) {
    if (k > 1) {
      d = z + d * rz1/rz2;
    }

    // w <- Jacobian()
    Jacobian(d, Omega12, P, w);

    double denom = arma::accu(d % w);
    double normr = arma::norm(rv, 2);

    if (denom <= 0) {
      p = d / arma::norm(d, 2);
      return;
    } else {
      double alpha = rz1 / denom;
      p += alpha * d;
      rv -= alpha * w;
    }

    normr = arma::norm(rv, 2);
    if (normr <= tolb) {
      return;
    } else {
      z = rv / c;
      rz2 = rz1;
      rz1 = arma::accu(rv % z);
    }
  }
  return;
}


void PCA(const arma::vec& lambda, const arma::mat& P, arma::mat& X) {

  int n = P.n_cols;
  int r = arma::accu(lambda > 0);
  int s = n - r;

  if (r == 0) {
    X.zeros();
  } else if (r == n) {
    // X = X
    return;
  } else if (r == 1) {
    X = lambda(0) * lambda(0) * (P.col(0) * P.col(0).t());
    return;
  } else if (r <= s) {
    arma::mat P1 = P.head_cols(r);
    arma::vec l1 = arma::sqrt(lambda.head(r));
    arma::mat P1l1 = P1.each_row() % l1.t();
    X = P1l1 * P1l1.t();
    return;
  } else { // (r > s)
    arma::mat P2 = P.tail_cols(s);
    arma::vec l2 = arma::sqrt(-1.0 * lambda.tail(s));
    arma::mat P2l2 = P2.each_row() % l2.t();
    X += P2l2 * P2l2.t();
    return;
  }
  return;
}


void cov2cor(arma::mat& X) {
  arma::mat D = arma::pinv(arma::diagmat(arma::sqrt(arma::diagvec(X))));
  X = D * X * D;
  for (unsigned int i = 0; i < X.n_cols; ++i) {
    X(i, i) = 1.0;
  }
  X = arma::symmatu(X);
  return;
}

//' Calculate the nearest positive semi-definite correlation matrix
//'
//' @param G the input correlation matrix
//' @param tau A user-dependent tuning parameter that determines the accuracy
//'     of the final correlation matrix. Smaller values generally mean faster
//'     convergence
//' @param iter_outer the max number of iterations in the outer loop
//' @param iter_inner the max number of iterations in the inner loop
//' @param maxit Maximum number of iterations in the pre_cg routine
//' @param err_tol the error tolerance for the stopping criteria
//' @param precg_err_tol the error tolerance in the pre-conjugate gradient method
//' @param newton_err_tol the error tolerance in Newton's method
//' @export
// [[Rcpp::export]]
arma::mat cor_nearPSD(
    arma::mat G,
    double tau=1e-5,
    int iter_outer=200,
    int iter_inner=20,
    int maxit=200,
    double err_tol=1e-6,
    double precg_err_tol=1e-2,
    double newton_err_tol=1e-4) {

  int n = G.n_cols;

  // Make G symmetric correlation matrix
  G = arma::symmatu(G);
  G.diag().ones();

  arma::vec b(n, arma::fill::ones);
  if (tau > 0) {
    b -= tau;
    G -= tau*arma::eye(n, n);
  }
  arma::vec b0(b);

  // Set X, compute eigen decomp, and sort values/vectors
  arma::vec y(n, arma::fill::zeros);
  arma::mat X = G + arma::diagmat(y);
  arma::vec lambda;
  arma::mat P;
  arma::eig_sym(lambda, P, X);
  lambda = arma::reverse(lambda);
  P = arma::fliplr(P);

  arma::vec Fy(n);
  double f0;
  gradient(y, lambda, P, b0, f0, Fy);  // modifies f0 and Fy
  double f = f0;
  b = b0 - Fy;

  arma::mat Omega12;
  omega_mat(lambda, Omega12);  // modifies Omega12
  arma::vec x0(y);

  PCA(lambda, P, X);
  double val_G    = 0.5 * std::pow(arma::norm(G, 2), 2);
  double val_dual = val_G - f0;
  double val_obj  = 0.5 * std::pow(arma::norm(X - G, 2), 2);
  double gap      = (val_obj - val_dual) / (1 + std::abs(val_dual) + std::abs(val_obj));

  double normb = arma::norm(b, 2);
  double normb0 = arma::norm(b0, 2) + 1;
  double delta_nb = normb / normb0;

  arma::vec c(n, arma::fill::ones);
  arma::vec d(n, arma::fill::zeros);

  int k = 0;
  while ( (gap > err_tol) && (delta_nb > err_tol) && (k < iter_outer) ) {
    precond_matrix(Omega12, P, c);                      // modifies c
    pre_cg(b, c, Omega12, P, precg_err_tol, maxit, d);  // modifies d

    double slope = arma::accu((Fy - b0) % d);

    y = x0 + d;
    X = G + arma::diagmat(y);
    arma::eig_sym(lambda, P, X);
    lambda = arma::reverse(lambda);
    P = arma::fliplr(P);
    gradient(y, lambda, P, b0, f, Fy);  // modifies f and Fy

    int k_inner = 0;
    while ( (k_inner <= iter_inner) &&
            (f > f0 + newton_err_tol * slope * std::pow(0.5, k_inner) + 1.0e-6) ) {
      ++k_inner;
      y = x0 + d * std::pow(0.5, k_inner);
      X = G + diagmat(y);
      arma::eig_sym(lambda, P, X);
      lambda = arma::reverse(lambda);
      P = arma::fliplr(P);
      gradient(y, lambda, P, b0, f, Fy);  // modifies f0 and Fy
    }

    x0 = y;
    f0 = f;

    PCA(lambda, P, X);
    val_dual = val_G - f0;
    val_obj  = 0.5 * std::pow(arma::norm(X - G, 2), 2);
    gap      = (val_obj - val_dual) / (1 + std::abs(val_dual) + std::abs(val_obj));

    b = b0 - Fy;
    normb = arma::norm(b, 2);
    delta_nb = normb / normb0;

    omega_mat(lambda, Omega12);  // modifies Omega12

    ++k;
  }

  X += tau*arma::eye(n, n);
  cov2cor(X);
  return X;
}
