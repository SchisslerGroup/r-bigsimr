test_that("jax_rmvn generates random multivariate normal data", {

  skip_if_not(have_jax(), "jax is not available for testing")

  # Should fail for negative definite matrices
  nd_cor <- matrix(c(
    1.00, 0.82, 0.56, 0.44,
    0.82, 1.00, 0.28, 0.85,
    0.56, 0.28, 1.00, 0.22,
    0.44, 0.85, 0.22, 1.00
  ), 4, 4, byrow=TRUE)
  d <- ncol(nd_cor)

  # expect_error(jax_rmvn(10, rep(0, d), nd_cor))

  n       <- 1e6
  pd_cor  <- cor_nearPD(nd_cor)
  mu      <- stats::rnorm(d, 10, 10)
  x       <- jax_rmvn(n, mu, pd_cor)
  mu_hat  <- colMeans(x)
  cor_hat <- cor(x)

  for (i in 1:d) {
    expect_equal(mu[i], mu_hat[i], tolerance=0.01)
  }
  expect_equal(pd_cor, cor_hat, tolerance=0.01)
})
