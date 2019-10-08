#' Generate a random Pearson correlation matrix
#'
#' @param d A positive integer.
#' @param constant_rho A boolean value.
#' @return A \code{dxd} symmetric matrix with ones on the diagonal.
#' @examples
#' rcor(4)
#' rcor(10, TRUE)
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
#' @param shape The shape parameter of a negative binomial distribution.
#' @param id_margins A boolean value specifying if the margins should be identical.
#' @return A matrix with columns containing (in order) lambda, alpha, mean, and variance.
#' @examples
#' rnbinom_params(4)
#' rnbinom_params(10, shape = 40)
#' rnbinom_params(7, 44, TRUE)
rnbinom_params <- function(d, shape = 100, id_margins = FALSE) {
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


#' Transforms a [multivariate]normal vector to a different marginal via a
#' uniform intermediate transformation.
#'
#' @param x A normal random vector.
#' @param param A list ccontaining the marginal and its parameters.
normal2marginal <- function(x, param) {
  do.call(what = paste0("q", param[[1]]),
          args = c(list(p = pnorm(x)), param[-1]))
}
