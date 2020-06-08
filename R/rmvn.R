.rmvn_jax <- function(seed, mu, rho, n) {

  # Any RNG in jax requires a key. If none is supplied by the user, then it
  # should be generated at random.
  key = jax$random$PRNGKey(sample(.Random.seed, 1))

  d <- ncol(rho)

  size     = reticulate::tuple(as.integer(n), as.integer(d))
  rho_np   = reticulate::np_array(rho)
  mu_np    = reticulate::np_array(mu)

  rho_chol = numpy$linalg$cholesky(rho_np)

  z      = numpy$random$normal(0.0, 1.0, size)

  z_d    = jax$device_put(z)
  rho_chol_d = jax$device_put(rho_chol)

  x_d = jax$numpy$matmul(z_d, rho_chol_d)

  x   = jax$device_get(x_d)

  reticulate::py_to_r(z)
}


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

