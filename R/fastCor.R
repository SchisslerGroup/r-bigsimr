#' Use the fastest methods available for computing each type of correlation
#'
#' @param x A matrix or vector
#' @param y A vector
#' @param method The type of correlation to compute
#' @export
fastCor <- function(x, y = NULL, method = c("pearson", "kendall", "spearman")) {

  method <- match.arg(method)

  stopifnot(method %in% c("pearson", "kendall", "spearman"))

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (method == "pearson") {
    if (is.null(y)) {
      coop::pcor(x)
    } else {
      coop::pcor(x, y)
    }

  } else if (method == "spearman") {
    if (is.null(y)) {
      coop::pcor(apply(x, 2, fastrank::fastrank_average))
    } else {
      coop::pcor(fastrank::fastrank_average(x),
                 fastrank::fastrank_average(y))
    }

  } else {
    if (is.null(y)) {
      pcaPP::cor.fk(x)
    } else {
      pcaPP::cor.fk(x, y)
    }

  }

}
