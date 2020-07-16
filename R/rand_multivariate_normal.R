#' Simulate random data from a multivariate normal distribution
#'
#' This function is designed to be 'bare bones' and does not check that the
#'   given covariance matrix is positive semi-definite.
#'
#' @param n A number of random vectors to be simulated.
#' @param mu A vector of length d, representing the mean.
#' @param sigma A symmetric positive definite covariance matrix (d x d)
#' @export
rmvn <- function(n, mu, sigma) {
  S    = jax$device_put(reticulate::np_array(sigma)) #
  m    = jax$device_put(reticulate::np_array(mu))
  key  = jax$random$PRNGKey(sample(.Random.seed, 1))
  size = reticulate::tuple(as.integer(n))
  z    = jax$random$multivariate_normal(key, m, S, size)
  reticulate::py_to_r(jax$device_get(z))
}

