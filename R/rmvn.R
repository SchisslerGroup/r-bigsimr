.rmvn_jax <- function(seed, mu, rho, n) {

  # Any RNG in jax requires a key. If none is supplied by the user, then it
  # should be generated at random.
  # if (is.null(seed)) {
  #   set.seed(seed)
  #   seed <- as.integer(sample(1:1000, 1))
  # }
  # key = jax$random$PRNGKey(as.integer(seed))
  #
  # rmvn <- jax$random$multivariate_normal

  size   = reticulate::tuple(n)

  rho_np = reticulate::np_array(rho)
  mu_np  = reticulate::np_array(mu)
  x      = numpy$random$multivariate_normal(mu_np, rho_np, size)

  # dev <- reticulate::py_to_r(jax$lib$xla_bridge$get_backend()$platform)
  #
  # mu    = jax$device_put(mu)
  # Sigma = jax$device_put(Sigma)
  #
  # if (dev == 'gpu') {
  #   x = rmvn(key, mu, Sigma, reticulate::tuple(list(n)))$block_until_ready()
  # } else {
  #   x = rmvn(key, mu, Sigma, reticulate::tuple(list(n)))
  # }
  #
  # x = jax$device_get(x)

  reticulate::py_to_r(x)
}
