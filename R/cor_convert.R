#' Convert from one correlation to another
#'
#' @param rho A correlation matrix
#' @param from Input correlation matrix
#' @param to Output correlation matrix
#' @export
cor_convert <- function(rho,
                        from = c("pearson", "spearman", "kendall"),
                        to = c("pearson", "spearman", "kendall")) {

  from <- match.arg(from)
  to   <- match.arg(to)

  # The integer mapping assigned here is arbitrary, but matches the ordering
  # found in `src/cor2cor.cpp`. It is important to either not change these
  # numbers or to ensure that they match in the other file.
  CASE <- switch (paste(from, to, sep = "_"),
                  "pearson_spearman" = 0L,
                  "pearson_kendall"  = 1L,
                  "spearman_pearson" = 2L,
                  "spearman_kendall" = 3L,
                  "kendall_pearson"  = 4L,
                  "kendall_spearman" = 5L,
                  -1L
  )

  if (is.matrix(rho)) {
    x <- .cor_convert_matrix(rho, CASE)
    x <- (x + t(x)) / 2
    diag(x) <- 1
    return(x)
  } else {
    return(.cor_convert_double(rho, CASE))
  }
}
