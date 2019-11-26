#' Generate a random Pearson correlation matrix
#'
#' @param d A positive integer.
#' @param constant_rho A boolean value.
#' @return A \code{dxd} symmetric matrix with ones on the diagonal.
#' @examples
#' rcor(4)
#' rcor(10, TRUE)
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


#' Generate negative binomial marginal parameters
#'
#' @param d A positive integer.
#' @param id_margins A boolean value specifying if the margins should be
#'   identical.
#' @return A matrix with columns containing (in order) lambda, alpha, mean,
#'   and variance.
#' @examples
#' rnbinom_params(4)
#' rnbinom_params(10)
#' rnbinom_params(7, TRUE)
#' @export
rnbinom_params <- function(d, id_margins = FALSE) {
  if (!id_margins) {
    probs <- runif(d, min = 0.3, max = 0.7)
    sizes <- sample(3:8, d, TRUE)
  } else {
    probs <- rep(runif(1, 0.3, 0.7), d)
    sizes <- rep(sample(3:8, 1, TRUE), d)
  }
  purrr::transpose(list(rep("nbinom", d), size = sizes, prob = probs))
}
