#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#define PERTURBATION 1.0e-9


// TODO: determine if precond_mat and pre_cg are working properly

void Jacobian(
        const arma::vec& x,
        const arma::mat& Omega12,
        const arma::mat& P,
        arma::vec& Ax)
{
    int n = P.n_cols;
    int r = Omega12.n_rows;
    int s = Omega12.n_cols;
    arma::mat P1 = P.head_cols(r);
    arma::mat P2 = P.tail_cols(s);

    Ax.zeros();
    arma::mat Omega(r, s);

    if (r == 0)
    {
        // TODO: this should probably throw an exception
    }
    else if (r == n)
    {
        Ax = x * (1.0 + PERTURBATION);
    }
    else if (r < s)
    {
        arma::mat H1 = P1 + x * arma::ones(1, r);
        arma::mat Omega = Omega12 % (H1.t() * P2);

        arma::mat HT(n, n, arma::fill::zeros);
        HT.head_cols(r) += P1 * P1.t() * H1 + P2 * Omega.t();
        HT.tail_cols(s) += P1 * Omega;

        Ax = arma::sum(P % HT, 1) + x * PERTURBATION;
    }
    else // (r >= s)
    {
        arma::mat H1 = P2 % (x * arma::ones(1, s));
        arma::mat Omega = 1.0 - Omega12 + P1.t() * H1;

        arma::mat HT(n, n, arma::fill::zeros);
        HT.head_cols(r) += P2 * Omega.t();
        HT.tail_cols(s) += P2 * H1.t() * P2 + P1 * Omega;

        Ax = arma::sum(P % HT, 1) + x * (1.0 + PERTURBATION);
    }
}


void gradient(
    const arma::vec& y,
    const arma::vec& lambda,
    const arma::mat& P,
    const arma::vec& b0,
    double& f,
    arma::vec& Fy)
{
    int n = P.n_cols;

    arma::vec lambda_nonneg(lambda);
    lambda_nonneg.elem(find(lambda < 0)).fill(0);

    arma::mat Q(P);
    for (int j = 0; j < n; j++)
    {
        Q.col(j) *= std::sqrt(lambda_nonneg(j));
    }

    for (int i = 0; i < n; i++)
    {
        Fy[i] = arma::accu(Q.row(i) % Q.row(i));
    }
    f = arma::accu( arma::pow(lambda_nonneg, 2) / 2.0 - arma::accu(b0 - y) );
}


void omega_mat(
    const arma::vec& lambda,
    arma::mat& Omega12)
{
    int n = lambda.n_elem;

    // number of strictly positive eigenvalues
    // lambda is assumed to be in descending order
    int r = 0;
    while(lambda(r) > 0 && r < n) r++;
    int s = n - r;

    if (r == 0)
    {
        Omega12.set_size(0,0);
    }
    else if (r == n)
    {
        Omega12.set_size(n, n);
        Omega12.ones();
    }
    else // 0 < r < n
    {
        Omega12.set_size(r, s);
        arma::vec lambda_s = lambda.tail(s);  // non-positive eigenvalues
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < s; j++)
            {
                Omega12(i, j) = lambda(i) / (lambda(i) - lambda_s(j));
            }
        }
    }
}


void precond_matrix(
    const arma::mat& Omega12,
    const arma::mat& P,
    arma::vec& c)
{

    int n = P.n_cols;
    int r = Omega12.n_rows;
    int s = Omega12.n_cols;

    c.ones();

    if (r == 0 || r == n)
    {
        // do not modify anything
        return;
    }

    arma::mat H = P % P;
    if (r < s)
    {
        arma::mat H12 = H.head_cols(r) * Omega12; // [n,s] = [n,r]*[r,s]
        c = arma::sum(H, 1); // c = rowsums(H)
        c %= c;              // c = c % c
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < s; j++)
            {
                c(i) += 2.0 * H12(i,j) * H(i,r+j);
            }
            c(i) = std::max(c(i), 1.0e-8);
        }
    }
    else
    {
        arma::mat H12 = H.tail_cols(s) * (1.0 - Omega12).t();   // [n,r] = [n,s]*[s,r]
        c = arma::sum(H.tail_cols(s), 1);
        c %= c;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < r; j++)
            {
                c(i) += 2.0 * H(i,j) * H12(i,j);
            }
            double alpha = accu(H.row(i));
            c(i) = -1.0 * c(i) + alpha*alpha;
            c(i) = std::max(c(i), 1.0e-8);
        }
    }
}


void pre_cg(
    const arma::vec& b,
    double tol,
    int maxit,
    const arma::vec& c,
    const arma::mat& Omega12,
    const arma::mat& P,
    arma::vec& p)
{
    int n = P.n_cols;
    int r = Omega12.n_rows;
    int s = Omega12.n_cols;

    arma::vec rv(b);
    double n2b = arma::norm(b, 2);

    double tolb = tol * n2b;
    p.zeros();
    arma::vec z = rv / c;
    double rz1 = arma::accu(rv % z);
    double rz2 = 1.0;
    arma::vec d(z);
    arma::vec w(n);
    for (int k = 1; k <= maxit; k++)
    {
        if (k > 1)
        {
            d = z + d * rz1/rz2;
        }

        Jacobian(d, Omega12, P, w); // modifies w

        double denom = arma::accu(d % w);
        double normr = arma::norm(rv, 2);

        if (denom <= 0)
        {
            p = d / arma::norm(d, 2);
            return;
        }
        else
        {
            double alpha = rz1 / denom;
            p += alpha * d;
            rv -= alpha * w;
        }

        z = rv / c;
        normr = arma::norm(rv, 2);
        if (normr <= tolb)
        {
            return;
        }
        else
        {
            rz2 = 0.0;
            rz1 = arma::accu(rv % z);
        }
    }
}


// [[Rcpp::export]]
arma::mat nearPDcor(arma::mat G)
{
    int n = G.n_cols;

    arma::vec b(n, arma::fill::ones);
    arma::vec b0(b);

    // Make G symmetric
    G = arma::symmatu(G);
    G.diag().ones();

    arma::vec y(n, arma::fill::zeros);
    arma::vec Fy(n, arma::fill::zeros);
    arma::vec x0(y);

    int k = 0;
    int iter_outer = 200;
    int iter_inner = 20;
    int maxit = 200;
    double tol = 1.0e-2;
    double err_tol = 1.0e-6;
    double sigma1 = 1.0e-4;

    arma::vec c(n, arma::fill::ones);
    arma::vec d(n, arma::fill::zeros);

    // Set X, compute eigen decomp, and sort values/vectors
    arma::mat X = G + arma::diagmat(y);
    arma::vec lambda;
    arma::mat P;
    if (!arma::eig_sym(lambda, P, X)) {
        std::cout << "First Eigen decomposition failed" << std::endl;
        return X;
    }
    lambda = arma::reverse(lambda);
    P = arma::fliplr(P);

    double f0;
    gradient(y, lambda, P, b0, f0, Fy);  // modifies f0 and Fy

    double f = f0;
    b = b0 - Fy;
    double normb = arma::norm(b, 2);

    arma::mat Omega12;
    omega_mat(lambda, Omega12);  // modifies Omega12

    x0 = y;

    while ( (normb > err_tol) && (k < iter_outer) )
    {
        precond_matrix(Omega12, P, c);  // modifies c
        pre_cg(b, tol, maxit, c, Omega12, P, d);  // modifies d

        double slope = arma::accu( (Fy - b0) % d );

        y = x0 + d;
        X = G + arma::diagmat(y);
        // arma::eig_sym(lambda, P, X);
        if (!arma::eig_sym(lambda, P, X)) {
            std::cout <<
                "Eigen decomp failed at k = " << k << std::endl;
            std::cout << "y = " << std::endl << y << std::endl;
            return X;
        }
        lambda = arma::reverse(lambda);
        P = arma::fliplr(P);
        gradient(y, lambda, P, b0, f0, Fy);  // modifies f0 and Fy

        int k_inner = 1;
        while ( (k_inner <= iter_inner) &&
                (f > f0 + sigma1 * std::pow(0.5, k_inner) * slope + 1.0e-5) )
        {
            y = x0 + std::pow(0.5, k_inner) * d;
            X = G + diagmat(y);
            // arma::eig_sym(lambda, P, X);
            if (!arma::eig_sym(lambda, P, X)) {
                std::cout <<
                    "Eigen decomp failed at k = " << k <<
                    " and k_inner = " << k_inner << std::endl;
                std::cout << "y = " << std::endl << y << std::endl;
                return X;
            }
            lambda = arma::reverse(lambda);
            P = arma::fliplr(P);
            gradient(y, lambda, P, b0, f0, Fy);  // modifies f0 and Fy
            k_inner++;
        }

        x0 = y;
        f0 = f;
        b = b0 - Fy;
        normb = arma::norm(b, 2);

        omega_mat(lambda, Omega12);  // modifies Omega12

        k++;
    }

    // number of strictly positive eigenvalues
    // lambda is assumed to be in descending order
    int r = 0;
    while(lambda(r) > 0 && r < n) r++;
    int s = n - r;

    if (r == 0)
    {
        X.zeros();
    }
    else if (r == n)
    {
        // do nothing
    }
    else if (r == 1)
    {
        X = lambda(0) * lambda(0) * (P.col(0) * P.col(0).t());
    }
    else if (r <= s)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < r; j++) {
                // all rows of column (r) are multiplied by lambda[r]
                P(i,j) *= std::sqrt(lambda(j));
            }
        }
        X = P.head_cols(r) * P.head_cols(r).t();
    }
    else // (r > s)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < s; j++)
            {
                // all rows of column (r+j) are multiplied by lambda[r+j]
                // I.e. last columns are multiplied by last lambdas
                P(i,r+j) *= std::sqrt(-1.0 * lambda(r+j));
            }
        }
        X = P.tail_cols(s) * P.tail_cols(s).t();
    }

    return X;
}
