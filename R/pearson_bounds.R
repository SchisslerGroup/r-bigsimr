#' Compute the theoretical lower and upper Pearson correlation bounds
#'
#' @param dA The first marginal distribution from the `distr6` package
#' @param dB The second marginal distribution from the `distr6` package
#' @param n The degree of the Hermite polynomial
#' @export
pearson_bounds <- function(dA, dB, n=7) {
  k <- 0:n
  a <- pm_getcoef(dA, n)
  b <- pm_getcoef(dB, n)

  mA <- mean(dA)
  mB <- mean(dB)
  sA <- distr6::stdev(dA)
  sB <- distr6::stdev(dB)

  c1 <- -mA * mB
  c2 <- 1 / (sA * sB)
  kab <- factorial(k) * a * b

  lower <- c1 * c2 + c2 * sum((-1)^k * kab)
  upper <- c1 * c2 + c2 * sum(kab)

  list(
    lower = clampcor(lower),
    upper = clampcor(upper)
  )

}
