#' Converts Pearson correlation to Spearman. Used for normal random variables.
#'
#' @param rho A a square symmetric Pearson correlation matrix.
#' @return A Spearman correlation matrix.
#' @export
convertPearson2Spearman <- function(rho) {
  tmp <- (6 / pi) * asin(rho / 2)
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(rho)
  return(tmp)
}


#' Converts Spearman correlation to Pearson. Used for normal random variables.
#'
#' @param rho A a square symmetric Spearman correlation matrix.
#' @return A Pearson correlation matrix.
#' @export
convertSpearman2Pearson <- function(rho) {
  tmp <- 2* sin( rho * (pi / 6))
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(rho)
  return(tmp)
}


#' Converts Pearson correlation to Kendall. Used for normal random variables.
#'
#' @param rho A a square symmetric Pearson correlation matrix.
#' @return A Kendall correlation matrix.
#' @export
convertPearson2Kendall <- function(rho) {
  tmp <- (2 / pi) * asin(rho)
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(rho)
  return(tmp)
}


#' Converts Kendall correlation to Pearson. Used for normal random variables.
#'
#' @param rho A a square symmetric Kendall correlation matrix.
#' @return A Pearson correlation matrix.
#' @export
convertKendall2Pearson <- function(rho) {
  tmp <- sin( rho * (pi / 2) )
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(rho)
  return(tmp)
}
