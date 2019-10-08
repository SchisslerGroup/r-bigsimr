#' Generate a random Pearson correlation matrix
#' 
#' @param d A positive integer.
#' @param constant_rho A boolean value.
#' @return A \code{dxd} symmetric matrix with ones on the diagonal.
#' @examples
#' rcorr(4)
#' rcorr(10, TRUE)
rcorr <- function(d, constant_rho = FALSE) {
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
#' @param shape The shape parameter of a negative binomial distribution.
#' @param id_margins A boolean value specifying if the margins should be identical.
#' @return A matrix with columns containing (in order) lambda, alpha, mean, and variance.
#' @examples
#' rnegbin_params(4)
#' rnegbin_params(10, shape = 40)
#' rnegbin_params(7, 44, TRUE)
rnegbin_params <- function(d, shape = 100, id_margins = FALSE) {
  if (!id_margins) {
    p <- runif(d, min = 0, max = 0.1)
    lambda <- (1 - p) / p
    alpha <- rgamma(d, shape = shape)
    nb_mean <- lambda * alpha
    nb_var <- ((1 - p) / p^2) * alpha
  } else {
    p <- rep(runif(1, min = 0, max = 0.1), d)
    lambda <- (1 - p) / p
    alpha <- rep(rgamma(1, shape = shape), d)
    nb_mean <- lambda * alpha
    nb_var <- ((1 - p) / p^2) * alpha
  }
  to_return <- rbind(lambda, alpha, nb_mean, nb_var)
  colnames(to_return) <- paste0("Var", 1:d)
  return(to_return)
}
