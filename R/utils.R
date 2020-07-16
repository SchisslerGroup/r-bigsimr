# find out whether the user has conda installed and visible
#' @importFrom reticulate conda_binary
have_conda <- function() {
  conda_bin <- tryCatch(
    reticulate::conda_binary("auto"),
    error = function(e)
      NULL
  )
  ! is.null(conda_bin)
}

#' @importFrom reticulate py_available
have_python <- function() {
  tryCatch(
    reticulate::py_available(initialize = TRUE),
    error = function(e)
      FALSE
  )
}

#' @importFrom reticulate py_module_available
have_numpy <- function() {
  reticulate::py_module_available("numpy")
}

#' @importFrom reticulate py_module_available
have_jax <- function() {
  reticulate::py_module_available("jax")
}
