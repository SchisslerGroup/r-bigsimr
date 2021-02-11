#' Pearson matching for the continuous-continuous case
#'
#' @param r The target correlation
#' @param dA The first marginal distribution
#' @param dB The second marginal distribution
#' @param n The degree of the Hermite polynomial
pm_continuous_case <- function(r, dA, dB, n) {
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

  coef <- numeric(n+1)
  for (k in 1:n) {
    coef[k+1] <- c2 * a[k+1] * b[k+1] * factorial(k)
  }
  coef[1] <- c1 * c2 + c2 * a[1] * b[1] - r

  rts <- solve_poly_pm_one(coef, r)
  if (!is.na(rts)) {
    return (rts)
  }

  warning("The target correlation is not feasible. Returning the match to the nearest bound instead.")
  lower <- c1 * c2 + c2 * sum((-1)^k * kab)
  upper <- c1 * c2 + c2 * sum(kab)

  if (r > 0) {
    pm_continuous_case(upper - 0.001, dA, dB, n)
  } else {
    pm_continuous_case(lower + 0.001, dA, dB, n)
  }

}

#' Pearson matching for the discrete-discrete case
#'
#' @param r The target correlation
#' @param dA The first marginal distribution
#' @param dB The second marginal distribution
#' @param n The degree of the Hermite polynomial
pm_discrete_case <- function(r, dA, dB, n) {
  sA <- distr6::stdev(dA)
  sB <- distr6::stdev(dB)
  minA <- distr6::inf(dA)
  minB <- distr6::inf(dB)
  maxA <- distr6::sup(dA)
  maxB <- distr6::sup(dB)
  maxA <- ifelse(is.infinite(maxA), dA$quantile(0.9999), maxA)
  maxB <- ifelse(is.infinite(maxB), dB$quantile(0.9999), maxB)

  # Support sets
  A <- seq(minA, maxA, by = 1.0)
  B <- seq(minB, maxB, by = 1.0)

  a <- c(-Inf, qnorm(dA$cdf(A))) # ||a|| = ||A|| + 1
  b <- c(-Inf, qnorm(dB$cdf(B))) # ||b|| = ||B|| + 1

  c2 <- 1 / (sA * sB)

  coef <- numeric(n+1)
  for (k in 1:n) {
    coef[k+1] <- Gn0d(A, B, a, b, c2, k) / factorial(k)
  }
  coef[1] <- -r

  rts <- solve_poly_pm_one(coef, r)
  if (!is.na(rts)) {
    return (rts)
  }

  warning("The target correlation is not feasible. Returning the match to the nearest bound instead.")
  bounds <- pearson_bounds(dA, dB, n)

  if (r > 0) {
    pm_discrete_case(bounds$upper - 0.001, dA, dB, n)
  } else {
    pm_discrete_case(bounds$lower + 0.001, dA, dB, n)
  }

}
