numpy <- NULL
jax <- NULL

.onLoad <- function(...) {
  # Check that a python environment exists and is available
  numpy <<- reticulate::import("numpy", delay_load = TRUE, convert = FALSE)
  jax <<- reticulate::import("jax", delay_load = TRUE, convert = FALSE)
}
