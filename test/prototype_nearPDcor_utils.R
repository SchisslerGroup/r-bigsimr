# Verified
set_omega_mat <- function(lambda) {
    n <- length(lambda)
    r <- sum(lambda > 0)
    s <- n - r

    if (r == 0) {
        return(matrix(0, 0, 0))
    } else if (r == n) {
        return(matrix(1, n, n))
    } else {
        M <- matrix(NA_real_, r, s)
        lambda_s <- tail(lambda, s)
        for (i in 1:r) {
            for (j in 1:s) {
                M[i,j] <- lambda[i] / (lambda[i] - lambda_s[j])
            }
        }
        return(M)
    }
}

# Verified
R_gradient <- function(y, lambda, P, b0) {
    n <- nrow(P)
    r <- sum(lambda > 0)
    if (r == 0) {
        Fy <- rep(0, n)
        f  <- 0
    } else {
        lambda <- ifelse(lambda < 0, 0, lambda)
        # column j is multiplied by lambda[j]
        Fy <- rowSums(sweep(P, 2, lambda, '*') * P)
        f  <- 0.5 * sum(lambda^2) - sum(b0 * y)
    }
    list(f = f, Fy = Fy)
}


# Modified
R_precond_matrix <- function(Omega12, P) {
    n <- nrow(P)
    r <- nrow(Omega12)
    s <- ncol(Omega12)

    cv <- rep(1.0, n)

    if (r == 0 || r == n) {
        return(cv)
    }

    H  <- t(P * P)
    # [r,n]
    H1 <- H[1:r, , drop=FALSE] # first r rows of H [r,n]
    H2 <- H[(r+1):n, , drop=FALSE]
    # crossprod(A,B) = A.T * B
    # tcrossprod(A,B) = A * B.T
    if (r < s) {
        H12 <- crossprod(H1, Omega12)    # [n,s] = [r,n].T*[r,s]
        cv  <- colSums(H1)^2 + 2*rowSums(H12 * t(H2))  # n = [r,n]->n + [n,s]%[s,n].T->n
    } else {
        H12 <- (1 - Omega12) %*% H2
        cv  <- colSums(H)^2 - colSums(H2)^2 - 2*rowSums(t(H1 * H12))
    }
    pmax(cv, 1e-8)
}


# Verified
R_pre_cg <- function(b, cv, Omega12, P, tol, maxit) {
    n <- nrow(P)
    r <- nrow(Omega12)
    s <- ncol(Omega12)

    n2b  <- norm(b, '2')
    tolb <- tol * n2b

    rv  <- b
    z   <- rv / cv
    rz1 <- sum(rv * z)
    rz2 <- 1.0

    p <- rep(0.0, n)
    d <- z
    z <- rv / cv

    for (k in 1:maxit) {
        if (k > 1) {
            d <- z + d * rz1/rz2;
        }

        w <- R_Jacobian(d, Omega12, P)

        denom <- sum(d * w)
        normr <- norm(rv, '2')
        if (denom <= 0) {
            return(d / norm(d, '2'))
        } else {
            alpha <- rz1 / denom
            p     <- p + alpha * d
            rv    <- rv - alpha * w
        }

        z <- rv / cv

        if (norm(rv, '2') <= tolb) {
            return(p)
        } else {
            rz2 <- rz1
            rz1 <- sum(rv * z)
        }
    }
    return(p)
}


# Assumed to be correct
R_Jacobian <- function(x, Omega12, P, PERTURBATION = 1e-9) {
    n <- nrow(P)
    r <- nrow(Omega12)
    s <- ncol(Omega12)

    if (0 < r && r < n) {
        P1 <- P[,     1:r, drop=FALSE]
        P2 <- P[, (r+1):n, drop=FALSE]
    }

    if (r == 0) {
        Ax <- rep(0, n)
    } else if (r == n) {
        Ax <- x * (1 + PERTURBATION)
    } else if (r < s) {
        H1    <- P1 + x %*% matrix(1, 1, r)
        Omega <- Omega12 * crossprod(H1, P2);

        HT1 <- tcrossprod(P1) %*% H1 + tcrossprod(P2, Omega)
        HT2 <- P1 %*% Omega

        Ax  <- rowSums(P * cbind(HT1, HT2)) + x * PERTURBATION
    } else {
        H1    <- P2 * (x %*% matrix(1, 1, s))
        Omega <- (1 - Omega12) + crossprod(P1, H1)

        HT1 <- P2 %*% t(Omega)
        HT2 <- tcrossprod(P2, H1) %*% P2 + (P1 %*% Omega)

        Ax  <- rowSums(P * cbind(HT1, HT2)) + x * (1 + PERTURBATION)
    }
    Ax
}

# This is the last step in the main.c program
R_PCA <- function(X, lambda, P) {
    n <- nrow(P)
    r <- sum(lambda > 0)
    s <- n - r

    if (r == 0) {
        X <- matrix(0, n, n)
    } else if (r == 1) {
        X <- lambda[1]^2 * tcrossprod(P[,1,drop=FALSE])
    } else if (r == n) {
        # X <- X
    } else if (r <= s) {
        P1   <- P[, 1:r, drop=FALSE]
        l1   <- sqrt(lambda[1:r])
        P1l1 <- sweep(P1, 2, l1, '*') # column j is multiplied by lambda
        X    <- tcrossprod(P1l1)
    } else {
        P2   <- P[, (r+1):n, drop=FALSE]
        l2   <- sqrt(-1*lambda[(r+1):n])
        P2l2 <- sweep(P2, 2, l2, '*')
        X    <- X + tcrossprod(P2l2)
    }

    X
}