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
