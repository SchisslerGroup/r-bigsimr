#' Compute a random positive definite correlation matrix
#'
#' @param d A positive integer number of dimensions
#' @param k A tuning parameter between 1 and d
#' @export
cor_randPD <- function(d, k=d) {
  cor_nearPD(cor_randPSD(d, k))
}


#' Compute a random positive semidefinite correlation matrix
#'
#' This method is generally faster than \code{\link{cor_randPD}}
#'
#' @param d A positive integer number of dimensions
#' @param k A tuning parameter between 1 and d
#' @export
cor_randPSD <- function(d, k=d) {
  stopifnot(d >= 1)
  stopifnot(1 <= k && k <= d)

  .cor_randPSD(d, k)
}
