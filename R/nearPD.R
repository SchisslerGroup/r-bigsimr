nearPD <- function (x,
                    eig.tol = 1e-06,
                    conv.tol = 1e-07,
                    posd.tol = 1e-08,
                    maxit = 100,
                    conv.norm.type = "I") {
  n <- ncol(x)

  D_S <- x
  D_S[] <- 0

  X <- x
  iter <- 0
  converged <- FALSE
  conv <- Inf
  while (iter < maxit && !converged) {
    Y <- X
    R <- Y - D_S
    e <- eigen(R, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    p <- d > eig.tol * d[1]

    if (!any(p))
      stop("Matrix seems negative semi-definite")

    Q <- Q[, p, drop = FALSE]
    X <- tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)

    D_S <- X - R

    diag(X) <- 1

    conv <- norm(Y - X, conv.norm.type) / norm(Y, conv.norm.type)
    iter <- iter + 1

    converged <- (conv <= conv.tol)
  }

  if (!converged)
    warning(gettextf("'nearPD()' did not converge in %d iterations", iter),
            domain = NA)

  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
    d[d < Eps] <- Eps
    Q <- e$vectors
    o.diag <- diag(X)
    X <- Q %*% (d * t(Q))
    D <- sqrt(pmax(Eps, o.diag) / diag(X))
    X[] <- D * X * rep(D, each = n)
  }

  diag(X) <- 1

  X

}
