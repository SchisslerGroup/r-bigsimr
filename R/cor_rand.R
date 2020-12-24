.rjm <- function(A, a) {
  b     <- dim(A)[1]
  idx   <- 2:(b-1)
  p1    <- A[idx,1,drop=FALSE]
  p3    <- A[idx,b,drop=FALSE]
  R2    <- A[idx,idx,drop=FALSE]
  Ri    <- solve(R2)
  rcond <- 2 * stats::rbeta(1, a, a) - 1
  t13   <- crossprod(p1, Ri) %*% p3
  t11   <- crossprod(p1, Ri) %*% p1
  t33   <- crossprod(p3, Ri) %*% p3
  c(t13 + rcond * sqrt((1 - t11) * (1 - t33)))
}


#' Compute a random positive definite correlation matrix
#'
#' @param d A positive integer number of dimensions
#' @param a A tuning parameter
#' @export
cor_randPD <- function(d, a=1.0) {
  if (d == 1) {
    return(matrix(1, 1, 1))
  } else if (d == 2) {
    p <- stats::runif(1)
    return(matrix(c(1, p, p, 1), 2, 2))
  } else {
    R <- diag(1, d, d)
    for (i in 1:(d-1)) {
      a0 <- a + (d - 2) / 2
      R[i, i+1] <- 2 * stats::rbeta(1, a0, a0) - 1
      R[i+1, i] <- R[i,i+1]
    }
    for (m in 2:(d-1)) {
      for (j in 1:(d-m)) {
        r_sub <- R[j:(j+m), j:(j+m),drop=FALSE]
        a0    <- a + (d - m - 1) / 2
        R[j, j+m] <- .rjm(r_sub, a0)
        R[j+m, j] <- R[j, j+m]
      }
    }
    R[] <- (R + t(R)) / 2
    return(R)
  }
}


#' Compute a random positive semidefinite correlation matrix
#'
#' This method is generally faster than \code{\link{cor_randPD}}
#'
#' @param d A positive integer number of dimensions
#' @param k A tuning parameter between 1 and d
#' @export
cor_randPSD <- function(d, k=d) {
  if (d == 1) {
    return(matrix(1, 1, 1))
  }
  stopifnot(1 <= k && k <= d)
  W  <- matrix(stats::rnorm(d * k), d, k)
  S  <- tcrossprod(W, W) + diag(stats::runif(d))
  S2 <- diag(1 / sqrt(diag(S)))
  R <- stats::cov2cor(S2 %*% S %*% S2)
  (R + t(R)) / 2
}
