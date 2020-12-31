#' Nearest positive (semi) definite correlation matrix
#'
#' @param G the input correlation matrix
#' @param tau smallest allowed eigenvalue. Set to 0 for semidefinite matrix
#' @param tol relative gap between the objective function value of the
#'     obtained solution and the objective function value of the global
#'     solution.
#' @export
cor_nearPD <- function(G, tau = 1e-6, tol = 1e-6) {
  n <- NCOL(G)
  if (n == 1) return(matrix(1,1,1))

  # Make G symmetric
  G[] <- Matrix::forceSymmetric(G, 'U')
  diag(G) <- 1.0

  # target diagonal (b = 1 for correlation matrix)
  b <- rep(1, n)
  if (tau > 0) {
    b <- b - tau
    diag(G) <- diag(G) - tau
  }
  b0 <- as.matrix(b)

  # initial dual value
  y  <- array(0, c(n, 1))
  Fy <- array(0, c(n, 1))

  # Maximum iterations
  iter_outer <- 200
  iter_inner <- 20
  iter_cg    <- 200

  # error tolerance in the gap
  error_tol <- max(1e-12, tol)
  # relative accuracy of the conjugate gradient method for solving the
  # Newton direction
  tol_cg <- 1e-2
  # tolerance in the line search of the Newton method
  tol_ls <- 1e-4

  x0 <- y
  cv <- array(1, c(n, 1))
  d  <- array(0, c(n, 1))

  X      <- G + diag(c(y))
  eig    <- eigen(X, symmetric = TRUE)
  P      <- eig$vectors
  lambda <- eig$values

  gradient <- npd_gradient(y, lambda, P, b0, n)
  f0       <- gradient$f
  Fy       <- gradient$Fy
  X        <- npd_pca(X, lambda, P, b0, n)

  val_G    <- 0.5 * norm(G, 'F') ^ 2
  val_dual <- val_G - f0
  val_obj  <- 0.5 * norm(X - G, 'F') ^ 2
  gap      <- (val_obj - val_dual) / (1 + abs(val_dual) + abs(val_obj))

  f <- f0
  b <- b0 - Fy

  norm_b     <- norm(b, '2')
  norm_b0    <- norm(b0, '2') + 1
  norm_b_rel <- norm_b / norm_b0

  Omega12 <- npd_omega_matrix(lambda, n)

  k <- 0
  while (gap > error_tol &&
         norm_b_rel > error_tol &&
         k < iter_outer) {

    cv <- npd_precondition_matrix(Omega12, P, n)
    d <- npd_pre_conjugate_gradient(b, tol_cg, iter_cg, cv, Omega12, P, n)

    slope <- sum((Fy - b0) * d)

    y <- x0 + d
    X      <- G + diag(c(y))
    eig    <- eigen(X, symmetric = TRUE)
    P      <- eig$vectors
    lambda <- eig$values

    gradient <- npd_gradient(y, lambda, P, b0, n)
    f <- gradient$f
    Fy <- gradient$Fy

    k_inner <- 0
    while (k_inner <= iter_inner &&
           f > f0 + tol_ls * 0.5 ^ k_inner * slope + 1e-6) {
      k_inner <- k_inner + 1

      y      <- x0 + 0.5^k_inner * d
      X      <- G + diag(c(y))
      eig    <- eigen(X, symmetric = TRUE)
      P      <- eig$vectors
      lambda <- eig$values

      gradient <- npd_gradient(y, lambda, P, b0, n)
      f0       <- gradient$f
      Fy       <- gradient$Fy
    }

    x0 <- y
    f0 <- f
    X  <- npd_pca(X, lambda, P, b0, n)
    val_dual <- val_G  - f0
    val_obj  <- 0.5 * norm(X - G, 'F')^2
    gap      <- (val_obj - val_dual) / (1 + abs(val_dual) + abs(val_obj))

    b <- b0 - Fy
    norm_b <- norm(b, '2')
    norm_b_rel <- norm_b / norm_b0

    Omega12 <- npd_omega_matrix(lambda, n)

    k <- k + 1
  }

  X[] <- Matrix::forceSymmetric(X, 'U')
  diag(X) <- 1.0
  X
}


npd_gradient <- function(y, lambda, P, b, n) {
  r <- sum(lambda > 0)

  if (r == 0) {
    Fy = matrix(0, n, 1)
    f = 0
  } else {
    P1 <- P[, 1:r, drop=FALSE]
    l1 <- lambda[1:r]

    Fy <- matrix(rowSums(sweep(P1, 2, l1, '*') * P1), n, 1)
    f  <- 0.5 * sum(l1 * l1) - sum(b * y)
  }

  list(f = f, Fy = Fy)
}


npd_pca <- function(X, lambda, P, b, n) {
  r <- sum(lambda > 0)
  s <- n - r

  if (n == 1) {
    return(b)
  }

  if (r == 0) {
    X <- mat.or.vec(n, n)
  } else if (r == n) {
    # X <- X
  } else if (r == 1) {
    X <- lambda[1] * tcrossprod(P[,1])
  } else if (r == 2) {
    P1   <- P[,1:r]
    l1   <- sqrt(lambda[1:r])
    P1l1 <- sweep(P1, 2, l1, '*')
    X    <- tcrossprod(P1l1)
  } else { # (r > s)
    P2   <- P[,(r+1):n, drop=FALSE]
    l2   <- sqrt(-lambda[(r+1):n])
    P2l2 <- sweep(P2, 2, l2, '*')
    X <- X + tcrossprod(P2l2)
  }

  # make the diagonal vector to be b without changing PSD property
  d <- diag(X)
  d <- pmax(d, b)
  X <- X - diag(diag(X)) + diag(d)
  d <- (b / d)^0.5
  d2 <- tcrossprod(d)
  X <- X * d2

  X
}


npd_omega_matrix <- function(lambda, n) {
  r <- sum(lambda > 0)
  s <- n - r

  if (r == 0) {
    return (matrix(0, 0, 0))
  } else if (r == n) {
    return (matrix(1, n, n))
  } else {
    l1 <- lambda[1:r]
    l2 <- lambda[(r+1):n]

    Omega12 <- matrix(1, r, s)
    Omega12[] <- apply(Omega12, 2, function(x) x * l1)
    Omega12[] <- sweep(Omega12, 2, l2, function(x, y) x / (x - y))

    return (Omega12)
  }
}


npd_jacobian_matrix <- function(d, Omega12, P, n) {
  r <- NROW(Omega12)
  s <- n - r

  perturbation <- 1e-10

  if (r == 0) {
    return (matrix(0, n, 1))
  }

  if (r == n) {
    return ((1 + perturbation) * d)
  }

  P1 <- P[,1:r, drop=FALSE]
  P2 <- P[,(r+1):n, drop=FALSE]

  O12 <- Omega12 * (t(P1) %*% sweep(P2, 1, d, '*'))
  PO <- P1 %*% O12

  hh <- 2 * rowSums(PO * P2)

  if (r <= s) {
    PP1 <- tcrossprod(P1)
    tmp <- apply(PP1, 1, function(x) sum(x^2 * d))
    Vd  <- tmp + hh + perturbation * d
  } else { # r > s
    PP2 <- tcrossprod(P2)
    tmp <- apply(PP2, 1, function(x) sum((x^2) * d))
    Vd  <- d + tmp + hh - (2 * d * diag(PP2)) + perturbation * d
  }

  Vd
}


npd_precondition_matrix <- function(Omega12, P, n) {
  r <- NROW(Omega12)
  s <- n - r

  if (r == 0 || r == n) return (matrix(1, n, 1))

  H <- t(P)
  H <- H * H
  H1 <- H[1:r,,drop=FALSE]
  H2 <- H[(r+1):n,,drop=FALSE]

  if (r < s) {
    H12 <- crossprod(H1, Omega12)
    cv <- colSums(H1)^2 + 2 * rowSums(H12 * t(H2))
  } else {
    H12 <- (1 - Omega12) %*% H2
    cv <- colSums(H)^2 - colSums(H2)^2 - 2 * rowSums(t(H1 * H12))
  }

  cv[cv < 1e-8] <- 1e-8
  cv
}


npd_pre_conjugate_gradient <- function(b, tol, iter_cg, cv, Omega12, P, n) {
  n2b   <- norm(b, '2')
  tol_b <- tol * n2b

  p   <- matrix(0, n, 1)
  r   <- b

  z   <- r / cv
  d   <- z
  rz1 <- sum(r * z)
  rz2 <- 1

  for (k in 1:iter_cg) {
    if (k > 1) {
      d <- z + (rz1 / rz2) * d
    }

    w <- npd_jacobian_matrix(d, Omega12, P, n) # w=V(y)*d

    denom <- sum(d * w)

    if (denom <= 0) {
      return (d / norm(d))
    } else {
      alpha <- rz1 / denom
      p <- p + alpha * d
      r <- r - alpha * w
    }

    if (norm(r, '2') <= tol_b) {
      return (p)
    }

    z <- r / cv
    rz2 <- rz1
    rz1 <- sum(r * z)
  }
  p
}
