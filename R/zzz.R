jax <- NULL
numpy <- NULL

.onLoad <- function(...) {
  # Check that a python environment exists and is available
  jax <<- reticulate::import("jax",
                             delay_load = TRUE,
                             convert = FALSE)
  numpy <<- reticulate::import("numpy",
                               as = "numpy",
                               delay_load = TRUE,
                               convert = FALSE)

  if (.Platform$OS.type == "windows" || !have_jax()) {
    options(use_jax = FALSE)
  } else {
    options(use_jax = getOption("use_jax", TRUE))
  }

}
