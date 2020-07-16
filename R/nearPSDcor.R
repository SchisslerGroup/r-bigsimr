.npsd_gradient <- function(y, eigvals, eigvecs, b0, n) {
  r <- sum(eigvals > 0)
  if (r == 0) {
    list(f = 0, Fy = rep(0, n))
  } else {
    eigvals[eigvals < 0] <- 0
    Fy <- rowSums(sweep(eigvecs, 2, eigvals, "*") * eigvecs)
    f  <- 0.5 * sum(eigvals^2) - sum(b0 * y)
    list(f = f, Fy = Fy)
  }
}

.npsd_pca <- function(X, eigvals, eigvecs, n) {
  r <- sum(eigvals > 0)
  s <- n - r
  if (r == 0) {
    matrix(0, n, n)
  } else if (r == 1) {
    eigvals[1]^2 * tcrossprod(eigvecs[,1,drop=FALSE], eigvecs[,1,drop=FALSE])
  } else if (r == n) {
    X
  } else if (r <= n) {
    P1   <- eigvecs[,1:r,drop=FALSE]
    l1   <- sqrt(utils::head(eigvals, r))
    P1l1 <- sweep(P1, 2, l1, "*")
    tcrossprod(P1l1, P1l1)
  } else {
    P2   <- eigvecs[,(r+1):n,drop=FALSE]
    l2   <- sqrt(-utils::tail(eigvals, s))
    P2l2 <- sweep(P2, 2, l2, "*")
    X + tcrossprod(P2l2, P2l2)
  }
}

.npsd_pre_cg <- function(b, cv, Omega0, eigvecs, eps, N, n) {
  eps_b <- eps * norm(b, "2")

  rv  <- b
  zv  <- rv / cv
  rz1 <- sum(rv * zv)
  rz2 <- 1.0

  p <- rep(0, n)
  d <- zv

  for (k in 1:N) {
    if (k > 1) {
      d <- zv + d * rz1 / rz2
    }

    w <- .npsd_jacobian(d, Omega0, eigvecs, n)

    denom <- sum(d * w)
    normr <- norm(rv, "2")
    if (denom <= 0) {
      return(d / norm(d, "2"))
    } else {
      a  <- rz1 / denom
      p  <- p + a * d
      rv <- rv - a * w
    }

    if (norm(rv, "2") <= eps_b) {
      return(p)
    } else {
      zv <- rv / cv
      rz2 <- rz1
      rz1 <- sum(rv * zv)
    }
  }
  p
}

.npsd_precond_matrix <- function(Omega0, eigvecs, n) {
  r <- dim(Omega0)[1]
  s <- dim(Omega0)[2]
  # r + s == n

  cv <- rep(1.0, n)

  if (r == 0 || r == n) {
    return(cv)
  }

  H  <- t(eigvecs^2)
  H1 <- H[1:r,,drop=FALSE]
  H2 <- H[(r+1):n,,drop=FALSE]

  if (r < s) {
    H12 <- crossprod(H1, Omega0)
    cv  <- colSums(H1)^2 + 2 * rowSums(tcrossprod(H12, H2))
  } else {
    H12 <- (1.0 - Omega0) %*% H2
    cv  <- colSums(H)^2 - colSums(H2)^2 - 2 * colSums(H1 * H12)
  }
  cv[cv < 1e-8] <- 1e-8
  cv
}

.npsd_set_omega <- function(eigvals, n) {
  r <- sum(eigvals > 0)
  s <- n - r

  if (r == 0) {
    matrix(0, 0, 0)
  } else if (r == n) {
    matrix(1.0, n, n)
  } else {
    M  <- matrix(0, r, s)
    lr <- utils::head(eigvals, r)
    ls <- utils::tail(eigvals, s)
    for (i in 1:r) {
      for (j in 1:s) {
        M[i,j] <- lr[i] / (lr[i] - ls[j])
      }
    }
    M
  }
}

.npsd_jacobian <- function(x, Omega0, eigvecs, n, PERTURB=1e-9) {
  r <- dim(Omega0)[1]
  s <- dim(Omega0)[2]
  # r + s == n

  if (0 < r && r < n) {
    P1 <- eigvecs[, 1:r, drop=FALSE]
    P2 <- eigvecs[, (r+1):n, drop=FALSE]
  }

  if (r == 0) {
    rep(0, n)
  } else if (r == n) {
    x * (1 + PERTURB)
  } else if (r < s) {
    H1    <- diag(x) %*% P1
    Omega <- Omega0 * crossprod(H1, P2)

    HT1 <- tcrossprod(P1, P1) %*% H1 + tcrossprod(P2, Omega)
    HT2 <- P1 %*% Omega

    rowSums(eigvecs * cbind(HT1, HT2)) + x * PERTURB
  } else {
    H2    <- diag(x) %*% P2
    Omega <- (1 - Omega0) * crossprod(P1, H2)

    HT1 <- tcrossprod(P2, Omega)
    HT2 <- tcrossprod(P2, H2) %*% P2 + P1 %*% Omega

    x * (1 + PERTURB) - rowSums(eigvecs * cbind(HT1, HT2))
  }
}

#' Compute the nearest positive semidefinite correlation matrix
#'
#' @param R the input correlation matrix
#' @param tau some non-negative threshold for the smallest eigenvalue
#' @param iter_outer the max number of iterations in the outer loop
#' @param iter_inner the max number of iterations in the inner loop
#' @param N the max number of iterations in the pre-conjugate gradient method
#' @param err_tol the error tolerance for the stopping criteria
#' @param precg_err_tol the error tolerance in the pre-conjugate gradient method
#' @param newton_err_tol the error tolerance in Newton's method
#' @export
nearestPSDcor <- function(R, tau=1e-5, iter_outer=200L, iter_inner=20L, N=200L,
                          err_tol=1e-6, precg_err_tol=1e-2, newton_err_tol=1e-4) {
  n <- dim(R)[1]

  R <- (R + t(R)) / 2

  b <- rep(1, n)
  if (tau > 0) {
    b <- b - tau
    R <- R - diag(tau, n)
  }
  b0 <- b

  y <- rep(0, n)
  X <- R + diag(y)
  eig <- eigen(X)

  ret <- .npsd_gradient(y, eig$values, eig$vectors, b0, n)
  f0  <- ret$f
  Fy  <- ret$Fy
  f   <- f0
  b   <- b0 - Fy

  Omega0 <- .npsd_set_omega(eig$values, n)
  x0     <- y

  X        <- .npsd_pca(X, eig$values, eig$vectors, n)
  val_R    <- 0.5 * norm(R, "2")^2
  val_dual <- val_R - f0
  val_obj  <- 0.5 * norm(X - R, "2")^2
  gap      <- (val_obj - val_dual) / (1 + abs(val_dual) + abs(val_obj))

  normb  <- norm(b, "2")
  normb0 <- norm(b0, "2") + 1
  del_nb <- normb / normb0

  k <- 0L
  while ((gap > err_tol) && (del_nb > err_tol) && (k < iter_outer)) {
    cv <- .npsd_precond_matrix(Omega0, eig$vectors, n)
    d  <- .npsd_pre_cg(b, cv, Omega0, eig$vectors, precg_err_tol, N, n)

    slope <- sum((Fy - b0) * d)

    y   <- x0 + d
    X   <- R + diag(y)
    eig <- eigen(X)
    ret <- .npsd_gradient(y, eig$values, eig$vectors, b0, n)
    f   <- ret$f
    Fy  <- ret$Fy

    k_inner <- 0L
    while ((k_inner <= iter_inner) &&
           (f > f0 + newton_err_tol * slope * 0.5^k_inner + 1e-6)) {
      k_inner <- k_inner + 1L
      y <- x0 + d * 0.5^k_inner
      X   <- R + diag(y)
      eig <- eigen(X)
      ret <- .npsd_gradient(y, eig$values, eig$vectors, b0, n)
      f   <- ret$f
      Fy  <- ret$Fy
    }

    x0 <- y
    f0 <- f

    X        <- .npsd_pca(X, eig$values, eig$vectors, n)
    val_dual <- val_R - f0
    val_obj  <- 0.5 * norm(X - R, "2")^2
    gap      <- (val_obj - val_dual) / (1 + abs(val_dual) + abs(val_obj))
    b        <- b0 - Fy
    normb    <- norm(b, "2")
    del_nb   <- normb / normb0

    Omega0 <- .npsd_set_omega(eig$values, n)

    k <- k + 1L
  }

  stats::cov2cor(X + diag(tau, n))
}
