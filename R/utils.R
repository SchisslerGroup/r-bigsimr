#' A helper function for creating a list of marginal distributions
#'
#' Similar to `alist`, except that the right hand side of arguments are
#'     evaluated.
#'
#' @param ... A list of marginal quantile functions with their named arguments
#' @examples
#' x <- rnorm(500, 3.14, 1.2)
#' y <- rlnorm(500, 2.5, 1.0)
#' z <- rexp(500, 2.718)
#'
#' m <- mlist(
#'   qnorm(mean = mean(x), sd = sd(x)),
#'   qlnorm(meanlog = mean(log(y)), sdlog = sd(log(y))),
#'   qexp(rate = (1 / mean(z)))
#' )
#'
#' # will not work
#' mlist(a = 10, b = 3)
#' @export
mlist <- function(...) {
  args <- rlang::exprs(...)
  lapply(args, function(m) {
    rlang::call2(rlang::call_name(m), !!!lapply(rlang::call_args(m), eval))
  })
}


#' Check if conda is available
#' @importFrom reticulate conda_binary
#' @export
have_conda <- function() {
  conda_bin <- tryCatch(
    reticulate::conda_binary("auto"),
    error = function(e)
      NULL
  )
  ! is.null(conda_bin)
}

#' Check if python is available
#' @importFrom reticulate py_available
#' @export
have_python <- function() {
  tryCatch(
    reticulate::py_available(initialize = TRUE),
    error = function(e)
      FALSE
  )
}

#' Check if numpy is available
#' @importFrom reticulate py_module_available
#' @export
have_numpy <- function() {
  reticulate::py_module_available("numpy")
}

#' Check if jax is available
#' @importFrom reticulate py_module_available
#' @export
have_jax <- function() {
  reticulate::py_module_available("jax")
}
