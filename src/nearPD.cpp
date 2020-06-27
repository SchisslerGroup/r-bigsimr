#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#define PERTURBATION 1.0e-9

inline void Jacobian_matrix(
        const arma::vec& x,
        const arma::mat& Omega12,
        const arma::mat& P,
        arma::mat& Ax
        ) {

    int n = P.n_rows;
    int r = Omega12.n_rows;
    int s = Omega12.n_cols;

    arma::mat Omega(r, s);

    if (r == 0) {
    } else if (r == n) {
        Ax = x * (1.0 + PERTURBATION);
    } else if (r < n / 2.0) {
        arma::mat H1(n, r);
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < r; j++) {
                H1(i, j) = x[i] * P(i, j);
            }
        }

        // [r,s] = [r,s] % [r,n]*[n,s]
        Omega = Omega12 % (H1.t() * P.tail_cols(s));

        // P1 * P1.T * H1
        // [n,r] = [n,r]*[r,n]*[n,r]
        arma::mat HT = P.head_cols(r) * P.head_cols(r).t() * H1;

        // [n,r] = [n,s]*[s,r]
        HT += P.tail_cols(s) * Omega.t();
        HT.resize(n,n);

        // [n,s] = [n,r]*[r,s]
        HT.tail_cols(s) += P.head_cols(r) * Omega;

        for (int i = 0; i < n; i++) {
            Ax[i] += arma::sum(P.row(i) % HT.row(i));
            Ax[i] += x[i] * PERTURBATION;
        }

    } else {
        arma::mat H1(n, s);
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < s; j++) {
                H1(i, j) = x[i] * P(i, j+r);
            }
        }

        Omega = arma::ones<arma::mat>(r, s) - Omega12;

        // [r,s] = [r,n]*[n,s]
        Omega %= P.head_cols(r).t() * H1;

        arma::mat HT = arma::zeros<arma::mat>(n, n);

        // [n,r] = [n,s]*[s,r]
        HT.head_cols(r) += P.tail_cols(s) * Omega.t();

        // [n,s] = [n,s]*[s,n]*[n,s]
        HT.tail_cols(s) += P.tail_cols(s) * H1.t() * P.tail_cols(s);
        HT.tail_cols(s) += P.head_cols(r) * Omega;

        for (int i = 0; i < n; i++) {
            Ax[i] -= arma::sum(P.row(i) % HT.row(i));
            Ax[i] += x[i] * (1.0 + PERTURBATION);
        }
    }
}

inline void pre_cg(
    arma::vec b,
    double tol,
    int maxit,
    arma::vec c,
    arma::mat& Omega12,
    arma::mat& P,
    arma::vec& p,
    int* flag,
    double* relres,
    int* iterk) {

    int n = P.n_rows;

    arma::vec r(b);

    double n2b = arma::norm(b, 2);
    double tolb = tol * n2b;

    *flag   = 1;
    *iterk  = 0;
    *relres = 1000;

    arma::vec z = r / c;

    double rz1 = arma::sum(r % z); // sum of elementwise multiplication
    double rz2 = 1.0;

    arma::vec d(z);

    // Ax = w
    arma::vec w = arma::zeros<arma::vec>(n);

    for (int k = 1; k <= maxit; k++) {

        *iterk = k;

        if (k > 1) {
            d = z + d * rz1 / rz2;
        }

        Jacobian_matrix(d, Omega12, P, w);

        double denom = arma::sum(d % w);
        double normr = arma::norm(r, 2);

        *relres = normr / n2b;

        if (denom <= 0) {
            double normd = arma::norm(d, 2);
            p = d / normd;
            return;
        } else {
            double alpha = rz1 / denom;
            p += alpha * d;
            r -= alpha * w;
        }

        z = r / c;
        normr = arma::norm(r, 2);

        if (normr <= tolb) {
            *iterk = k;
            *relres = normr / n2b;
            *flag = 0;
            return;
        }

        rz2 = rz1;
        rz1 = arma::sum(r % z);
    }
    return;
}

inline void precond_matrix(arma::mat& Omega12, arma::mat& P, arma::vec& c) {
    int r = Omega12.n_rows;
    int s = Omega12.n_cols;
    int n = P.n_rows;

    c.ones();

    if (r == n || r == 0) return;

    arma::mat H = P % P;


    if (0 < r && r < n/2) {
        arma::mat H12 = H.head_rows(r).t() * Omega12;
        for (int i = 0; i < n; i++) {
            c[i] = 0;

            for (int j = 0; j < r; j++) {
                c[i] += H(i, j);
            }
            c[i] *= c[i];

            for (int j = 0; j < s; j++) {
                c[i] += 2 * H12(i, j) * H12(i, j+r);
            }

            c[i] = std::max(c[i], 1.0e-8);
        }
    } else { // r >= n/2
        arma::mat Omega = arma::ones(arma::size(Omega12)) - Omega12;
        arma::mat H12 = Omega * H.tail_rows(s);

        for (int i = 0; i < n; i++) {
            c[i] = 0;
            for (int j = 0; j < s; j++) {
                c[i] += H(i, j+r);
            }
            c[i] *= c[i];

            for (int j = 0; j < r; j++) {
                c[i] += 2.0 * H(i, j) * H12(i, j);
            }

            double alpha = arma::sum(H.row(i));
            c[i] = -1*c[i] + alpha * alpha;

            c[i] = std::max(c[i], 1.0e-8);
        }
    }
}

inline void gradient(
    arma::vec& y,
    arma::vec& lambda,
    arma::mat& P,
    arma::vec& b0,
    double& f,
    arma::vec& Fy) {

    f = 0;
    int n = P.n_rows;

    arma::vec tmp_lambda = lambda;
    tmp_lambda.elem(find(lambda < 0)).fill(0);
    arma::mat Q = P % arma::ones(n) * sqrt(tmp_lambda).t();

    for (int i = 0; i < n; i++)
        Fy[i] = arma::sum( Q.row(i) % Q.row(i) );

    f = arma::sum(tmp_lambda % tmp_lambda) / 2 - arma::sum(b0 % y);
}

inline void omega_mat(arma::vec& lambda, int n, arma::mat& Omega12) {
    int r = 0;
    while (lambda[r] > 0 && r < n) r++;

    if (r == 0) {
        Omega12.resize(0, 0);
    } else if (r == n) {
        Omega12 = arma::ones(n, n);
    } else {
        Omega12.resize(r, n-r);
        Omega12.zeros();
        for (int j = 0; j < n-r; j++) {
            for (int i = 0; i < r; i++) {
                Omega12(i, j) = lambda[i] / (lambda[i] - lambda[r+j]);
            }
        }
    }
}


// [[Rcpp::export]]
arma::mat nearPD(arma::mat G) {
    int i, j;
    int n = G.n_rows;
    arma::mat X(n, n);
    arma::vec y(n);

    double tau   = 0;
    arma::vec b  = arma::ones<arma::vec>(n);

    G = arma::symmatu(G);
    G.diag() -= tau;

    arma::vec b0 = b - tau;

    arma::vec Res_b = arma::zeros<arma::vec>(300);

    // Initial point
    y.zeros();
    arma::vec Fy = arma::zeros<arma::vec>(n);

    int k = 0;
    int f_eval = 0;

    int iter_outer = 200;
    int iter_inner = 20;
    int maxit = 200;
    int iterk = 0;
    double tol = 1.0e-2;

    double error_tol = 1.0e-6;
    double sigma_1 = 1.0e-4;

    arma::vec x0(y);
    arma::vec c = arma::ones<arma::vec>(n);
    arma::vec d = arma::zeros<arma::vec>(n);

    arma::vec val_G = arma::sum(arma::trimatu(G) % arma::trimatu(G) / 2.0);
    X = G + arma::diagmat(y);
    arma::mat P(X);

    // Compute the Eigen decomposition of P
    // Store the values in lambda and the vectors in P
    arma::vec lambda;
    arma::eig_sym(lambda, P, P);
    lambda = arma::reverse(lambda); // Need eigenvalues in descending order
    P = arma::fliplr(P); // eigenvector columns need to be reversed as well

    // Computes the gradient and stores the results in P and f0
    double f0;
    gradient(y, lambda, P, b0, f0, Fy);

    double f = f0;
    f_eval += 1;
    b = b0 - Fy;
    double normb = arma::norm(b, 2);

    arma::vec initial_f = val_G - f0;

    arma::mat Omega12;
    omega_mat(lambda, n, Omega12);
    x0 = y;

    while (normb > error_tol && k < iter_outer) {
        precond_matrix(Omega12, P, c);

        int flag;
        double relres;

        pre_cg(b, tol, maxit, c, Omega12, P, d, &flag, &relres, &iterk);

        double slope = arma::sum( (Fy - b0) % d );
        y = x0 + d;

        X = G + arma::diagmat(y);
        P = X;
        arma::eig_sym(lambda, P, P);
        lambda = arma::reverse(lambda);

        gradient(y, lambda, P, b0, f0, Fy);

        int k_inner = 0;
        while (k_inner <= iter_inner &&
               f > f0 + sigma_1 * std::pow(0.5, k_inner) * slope + 1e-5) {
            k_inner++;
            y = x0 + std::pow(0.5, k_inner) * d;
            X = G;
            X.diag() = G.diag() + y;
            P = X;
            arma::eig_sym(lambda, P, P);
            lambda = arma::reverse(lambda);
            gradient(y, lambda, P, b0, f0, Fy);
        }

        f_eval += k_inner + 1;

        x0 = y;
        f0 = f;

        b = b0 - Fy;

        normb = arma::norm(b, 2);

        Res_b[k] = normb;
        omega_mat(lambda, n, Omega12);

        k++;
    }

    int r = 0;
    while (lambda[r] > 0 && r < n) r++;

    if (r == n) {
    } else if (r == 0) {
        X.zeros();
    } else if (r == 1) {
        X = lambda[0]*lambda[0] * P.col(0) * P.col(0).t();
    } else if (r <= n / 2.0) {
        P.head_cols(r) %= arma::sqrt(arma::ones(n) * lambda.head(r).t());
        X = P.head_cols(r) * P.head_cols(r).t();
    } else {
        P.tail_cols(n-r) %= arma::sqrt(arma::ones(n) * lambda.tail(n-r).t());
        X = P.tail_cols(n-r) * P.tail_cols(n-r).t() + X;
    }

    arma::vec final_f = val_G - f;
    arma::vec val_obj = arma::sum(arma::pow(X - G, 2) / 2);
    X.diag() += tau;
    return X;
}
