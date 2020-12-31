#' Simulate random data from a multivariate normal distribution
#'
#' This function is designed to be 'bare bones' and does not check that the
#'   given covariance matrix is positive definite.
#'
#' @param n The number of random vectors to be simulated.
#' @param mu A vector of length d, representing the mean.
#' @param sigma A symmetric positive definite covariance matrix (d x d)
#' @export
jax_rmvn <- function(n, mu, sigma) {
  S    = jax$device_put(reticulate::np_array(sigma))
  m    = jax$device_put(reticulate::np_array(mu))
  key  = jax$random$PRNGKey(sample(.Random.seed, 1))
  size = reticulate::tuple(as.integer(n))
  z    = jax$random$multivariate_normal(key, m, S, size)
  reticulate::py_to_r(jax$device_get(z))
}


# Internal function for random multivariate uniform
.rmvuu <- function(n, rho) {
  R = jax$device_put(reticulate::np_array(rho))
  m = jax$numpy$zeros(reticulate::tuple(ncol(rho)))

  key  = jax$random$PRNGKey(sample(.Random.seed, 1))
  size = reticulate::tuple(as.integer(n))

  Z = jax$random$multivariate_normal(key, m, R, size)
  U = jax$scipy$stats$norm$cdf(Z)$block_until_ready()

  reticulate::py_to_r(jax$device_get(U))
}


#' Generate random data from a multivariate uniform distribution
#'
#' @param n number of random vectors to be simulated
#' @param rho correlation matrix (d x d)
#' @param min either a single number or a vector of length d
#' @param max either a single number or a vector of length d
#' @export
jax_rmvu <- function(n, rho, min = 0, max = 1) {
  d <- ncol(rho)

  rho <- stats::cov2cor(rho)

  if (length(min) > 1) {
    if (length(min) != d)
      stop("'min' must either be length 1 or the same length as the number of dimensions.")
    l <- matrix(1, n, 1) %*% min
  } else {
    l <- matrix(min, n, d)
  }

  if (length(max) > 1) {
    if (length(max) != d)
      stop("'max' must either be length 1 or the same length as the number of dimensions.")
    u <- matrix(1, n, 1) %*% max
  } else {
    u <- matrix(max, n, d)
  }

  stopifnot(all(max > min))

  .rmvuu(n, rho) * (u - l) + l
}
