#' Gauss-Hermite quadrature - "Probabilists"
#'
#' Computes the sample points and weights for Gauss-Hermite quadrature. These
#'     sample points and weights will correctly integrate polynomials of degree
#'     \eqn{2 * n - 1} or less over the interval \eqn{[-Inf, Inf]} with the weight
#'     function \eqn{f(x) = \exp(-x^2/2)}.
#'
#' @param n Number of sample points and weights. It must be >= 1.
#' @return A list of nodes and weights
#' @export
numpy_hermegauss <- function(n = 1L) {
  if (!is.numeric(n) || n < 1)
    stop("'n' must be a positive integer.")

  if (!is.integer(n))
    n <- as.integer(n)

  ret <- reticulate::py_to_r(numpy$polynomial$hermite_e$hermegauss(n))
  names(ret) <- c("nodes", "weights")
  ret
}
