#' Transforms a [multivariate]normal vector to a different marginal via a
#' uniform intermediate transformation.
#'
#' @param x A normal random vector.
#' @param param A list ccontaining the marginal and its parameters.
normal2marginal <- function(x, param) {
  do.call(what = paste0("q", param[[1]]),
          args = c(list(p = pnorm(x)), param[-1]))
}
