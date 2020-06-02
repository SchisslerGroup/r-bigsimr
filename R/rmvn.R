.rmvn_jax <- function(seed, mu, Sigma, n) {

  # Any RNG in jax requires a key. If none is supplied by the user, then it
  # should be generated at random.
  if (is.null(seed)) {
    set.seed(seed)
    seed <- as.integer(sample(1:1000, 1))
  }

  rmvn <- jax$random$multivariate_normal

  Sigma = numpy$array(Sigma, dtype=numpy$float32)
  mu  = numpy$array(mu, dtype=numpy$float32)

  key = jax$random$PRNGKey(as.integer(seed))

  dev <- reticulate::py_to_r(jax$lib$xla_bridge$get_backend()$platform)

  # mu  = jax$device_put(mu)
  # Sigma = jax$device_put(Sigma)
  x = rmvn(key, mu, Sigma, reticulate::tuple(list(n)))$block_until_ready()
  x = jax$device_get(x)

  reticulate::py_to_r(x)
}
