# Solve a polynomial on the interval [-1, 1]
solve_poly_pm_one <- function(coeffs, target) {
  xs <- polyroot(coeffs)
  xs <- Re(xs[abs(Im(xs)) < .Machine$double.eps]) # Approximately real solutions
  xs <- xs[-1 <= xs & xs <= 1]

  if (length(xs) == 0) { # no real roots found in [-1, 1]
    return(NA)
  } else if (length(xs) == 1) { # one real root found
    return (Re(xs))
  } else {
    warning("More than one root found in the interval [-1, 1]. Returning the closest match.")
    message("The set of roots are [", paste(xs, collapse = ", "), "]")

    idx <- which.min(abs(xs - target))

    return (Re(xs[idx]))
  }
}

# Compute the Hermite polynomial of degree n at x
pm_Hn <- function(x, n) {
  2^(-n/2) * fastGHQuad::evalHermitePoly(x / sqrt(2), n)
}

# Get the coefficients used in Pearson matching
pm_getcoef <- function(margin, n=7) {
  m <- n + 4

  gh <- fastGHQuad::gaussHermiteData(m)
  t <- gh$x*sqrt(2)
  w <- gh$w

  X <- margin$quantile(pnorm(t))

  # Hermite(t, k) for k in 0:n
  H <- vapply(0:n, function(k) pm_Hn(t, k), FUN.VALUE = numeric(m))

  # factorial(k) for k in 0:n
  f <- factorial(0:n)

  colSums(sweep(sweep(H, 1, w*X, "*"), 2, f, "/")) / sqrt(pi)
}
