#' Computes the Hermite polynomial
#' @param x a numeric vector giving the values at which the Hermite polynomial should be evaluated.
#' @param n an integer giving the degree of the Hermite polynomials
#' @export
hermite <- function(x, n) {
  x[] <- .Call(`_bigsimr_hermite`, x, n)
  x
}
