jax <- NULL
numpy <- NULL

.onLoad <- function(...) {
  # Check that a python environment exists and is available
  jax <<- reticulate::import("jax",
                             delay_load = TRUE,
                             convert = FALSE)
  numpy <<- reticulate::import("jax.numpy",
                               as = "numpy",
                               delay_load = TRUE,
                               convert = FALSE)
}
