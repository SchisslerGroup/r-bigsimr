#' Converts Pearson correlation to Spearman. Used for normal random variables.
#'
#' @param rho A a square symmetric Pearson correlation matrix.
#' @return A Spearman correlation matrix.
#' @export
convertPearsonSpearman <- function(rho) {
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
convertSpearmanPearson <- function(rho) {
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
convertPearsonKendall <- function(rho) {
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
convertKendallPearson <- function(rho) {
  tmp <- sin( rho * (pi / 2) )
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(rho)
  return(tmp)
}


#' Transforms a [multivariate]normal vector to a different marginal via a
#' uniform intermediate transformation.
#'
#' @param x A normal random vector.
#' @param param A list ccontaining the marginal and its parameters.
normal2marginal <- function(x, param) {
  do.call(what = paste0("q", param[[1]]),
          args = c(list(p = pnorm(x)), param[-1]))
}

