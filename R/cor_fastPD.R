#' Return a positive definite correlation matrix that is close to the input matrix
#'
#' @param R The input correlation matrix
#' @param tau A minimum eigenvalue. Set to zero for positive semidefinite.
#' @export
cor_fastPD <- function(R, tau=1e-6) {
  n <- ncol(R)
  tau <- max(1e-12, tau)

  R[] <- Matrix::forceSymmetric(R)
  diag(R) <- 1 - tau

  eig <- eigen(R, symmetric = TRUE)
  R[] <- .pca(R, eig$values, eig$vectors)

  diag(R) <- diag(R) + tau
  R[] <- cov2cor(R)
  R[] <- Matrix::forceSymmetric(R)
  R
}


.pca <- function(R, l, P) {
  n <- ncol(R)
  r <- sum(l > 0)
  s <- n - r

  if (r == 0) {
    return (matrix(0.0, n, n))
  } else if (r == n) {
    return (R)
  } else if (r == 1) {
    return ((l[1] * l[1]) * (tcrossprod(P[,1,drop=FALSE])))
  } else if (r <= s) {
    P1 <- P[,1:r,drop=FALSE]
    l1 <- sqrt(l[1:r])
    P1l1 <- eachrow(P1, l1, "*")
    return (tcrossprod(P1l1))
  } else {
    P2 <- P[,(r+1):n,drop=FALSE]
    l2 <- sqrt(-l[(r+1):n])
    P2l2 <- eachrow(P2, l2, "*")
    return (R + tcrossprod(P2l2))
  }

}
