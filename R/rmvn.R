#' Simulate random data from a multivariate normal distribution
#'
#' This function is designed to be 'bare bones' and does not check that the
#'   given covariance matrix is positive semi-definite.
#'
#' @param n number of random vectors to be simulated.
#' @param mu vector of length d, representing the mean.
#' @param sigma covariance matrix (d x d). Alternatively is can be the cholesky
#'   decomposition of the covariance. In that case isChol should be set to TRUE.
#' @param isChol boolean set to true if sigma is the cholesky decomposition of
#'   the covariance matrix.
#' @export
rmvn <- function(n, mu, sigma, isChol = FALSE) {

  S    = jax$device_put(reticulate::np_array(sigma)) #
  m    = jax$device_put(reticulate::np_array(mu))
  key  = jax$random$PRNGKey(sample(.Random.seed, 1))

  if (isChol) {
    size = reticulate::tuple(as.integer(n), as.integer(ncol(sigma)))
    z    = jax$numpy$add(jax$numpy$matmul(jax$random$normal(key, size), S), m)
  } else {
    size = reticulate::tuple(as.integer(n))
    z    = jax$random$multivariate_normal(key, m, S, size)
  }

  reticulate::py_to_r(jax$device_get(z))

}

