#' Setup Distributions.jl
#'
#' This function initializes the Distributions package
#'
#' @export
Distributions_setup <- function() {
  JuliaCall::julia_install_package_if_needed("Distributions")
  JuliaCall::julia_library("Distributions")
  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(Distributions)),\"!\"=>\"_bang\"))")
  functions <- c(functions, "rand")
  Distributions <- JuliaCall::julia_pkg_import("Distributions", functions)
  Distributions
}
