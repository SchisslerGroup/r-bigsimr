#' Check if a given matrix is a valid
#'
#' @export
is_valid_correlation <- function(rho) {
  all(diag(rho) == 1) && matrixcalc::is.positive.semi.definite(rho)
}
