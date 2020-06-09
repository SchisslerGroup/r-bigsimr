.rmvuu <- function(n, rho) {

  R    = reticulate::np_array(rho)
  R_d  = jax$device_put(R)
  m_d  = jax$numpy$zeros(reticulate::tuple(ncol(rho)))
  key  = jax$random$PRNGKey(sample(.Random.seed, 1))
  size = reticulate::tuple(as.integer(n))
  Z_d  = jax$random$multivariate_normal(key, m_d, R_d, size)
  U_d  = jax$scipy$stats$norm$cdf(Z_d)
  U    = jax$device_get(U_d)

  reticulate::py_to_r(U)
}


#' Generate random data from a multivariate uniform distribution
#'
#' @param n number of random vectors to be simulated
#' @param rho correlation matrix (d x d)
#' @param min Either a single number or a vector of length d
#' @param max Either a single number or a vector of length d
#'
rmvu <- function(n, rho, min = 0, max = 1) {
  d <- ncol(rho)

  rho <- cov2cor(rho)

  if (length(min) > 1) {
    if (length(min) != d)
      stop("'min' must either be length 1 or the same length as the number of dimensions.")
    l <- matrix(1, n, 1) %*% min
  } else {
    l <- matrix(min, n, d)
  }

  if (length(max) > 1) {
    if (length(max) != d)
      stop("'min' must either be length 1 or the same length as the number of dimensions.")
    u <- matrix(1, n, 1) %*% max
  } else {
    u <- matrix(max, n, d)
  }

  stopifnot(all(max > min))

  # .rmvuu(n, rho) * (u - l) + l
  list(l, u)
}
