.rmvn_jax <- function(seed, mu, rho, n) {

  # Any RNG in jax requires a key. If none is supplied by the user, then it
  # should be generated at random.
  if (is.null(seed)) {
    set.seed(seed)
    seed <- as.integer(sample(1:1000, 1))
  }
  key = jax$random$PRNGKey(as.integer(seed))

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
