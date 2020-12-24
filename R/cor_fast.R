#' Use the fastest methods available for computing each type of correlation
#'
#' @param x A matrix or vector
#' @param y A vector
#' @param method The type of correlation to compute
#' @export
cor_fast <- function(x, y = NULL, method = c("pearson", "kendall", "spearman")) {

  method <- match.arg(method)

  stopifnot(method %in% c("pearson", "kendall", "spearman"))

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  storage.mode(x) <- "numeric"
  if (!is.null(y))
    storage.mode(y) <- "numeric"

  if (method == "pearson") {
    if (is.null(y)) {
      coop::pcor(x)
    } else {
      coop::pcor(x, y)
    }

  } else if (method == "spearman") {
    if (is.null(y)) {
      coop::pcor(apply(x, 2, fastrank_num_avg))
    } else {
      coop::pcor(fastrank_num_avg(x),
                 fastrank_num_avg(y))
    }

  } else {
    if (is.null(y)) {
      pcaPP::cor.fk(x)
    } else {
      pcaPP::cor.fk(x, y)
    }

  }

}


fastrank_num_avg <- function(x) {
  .Call("fastrank_num_avg_", x, PACKAGE = "bigsimr")
}
