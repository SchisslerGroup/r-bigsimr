#' Generate a random Pearson correlation matrix
#'
#' @param d A positive integer
#' @param constant_rho A boolean value
#' @return A \code{dxd} symmetric matrix with ones on the diagonal
#' @export
rcor <- function(d, constant_rho = FALSE) {
  if (constant_rho) {
    tmp_rho <- runif(n = 1, min = 0.25, max = 0.9)
    sigma <- matrix(data = tmp_rho, nrow = d, ncol = d)
    diag(sigma) <- 1
  } else {
    tmp <- matrix(rnorm(d^2), d, d)
    mcov <- tcrossprod(tmp, tmp)
    sigma <- cov2cor(mcov)
  }
  rownames(sigma) <- colnames(sigma) <- paste0("Var", 1:d)
  return(sigma)
}
