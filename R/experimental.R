#' Convert from one correlation to another
#' @param rho A correlation matrix
#' @param from Input correlation matrix
#' @param to Output correlation matrix
#' @export
cor2cor <- function(rho,
                    from = c("pearson", "spearman", "kendall"),
                    to = c("pearson", "spearman", "kendall")) {

  from <- match.arg(from)
  to   <- match.arg(to)

  CASE <- switch (paste(from, to, sep = "_"),
               "pearson_spearman" = 0L,
               "pearson_kendall"  = 1L,
               "spearman_pearson" = 2L,
               "spearman_kendall" = 3L,
               "kendall_pearson"  = 4L,
               "kendall_spearman" = 5L,
               -1L
  )

  CXX_cor2cor(rho, CASE)
}
