#' @importFrom reticulate py_to_r
.rmvn_jax <- function(seed, mu, Sigma, n) {

  # Any RNG in jax requires a key. If none is supplied by the user, then it
  # should be generated at random.
  if (is.null(seed)) {
    set.seed(seed)
    seed <- as.integer(sample(1:1000, 1))
  }

  rmvn <- jax$random$multivariate_normal

  Sigma = numpy$array(Sigma, dtype=numpy$float32)
  Sigma = jax$device_put(Sigma)

  mu  = numpy$array(mu, dtype=numpy$float32)
  mu  = jax$device_put(mu)

  key = jax$random$PRNGKey(as.integer(seed))
  x = rmvn(key, mu, Sigma, list(n))$block_until_ready()
  x = jax$device_get(x)
  reticulate::py_to_r(x)
}
