#' Check if a given correlation matrix is a valid
#'
#' A correlation matrix is valid if:
#' * The diagonal elements are all equal to 1
#' * It is symmetric
#' * It is positive semi-definite
#'
#' @param rho Input correlation matrix
#' @export
is_valid_correlation <- function(rho) {
  all(diag(rho) == 1) &&
    matrixcalc::is.positive.semi.definite(rho) &&
    matrixcalc::is.symmetric.matrix(rho)
}
