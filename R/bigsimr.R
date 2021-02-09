#' Setup bigsimr
#'
#' This function initializes Julia and the MvSim.jl package.
#' The first time will be long since it includes precompilation.
#' Additionally, this will install Julia and the required packages
#' if they are missing.
#'
#' @param pkg_check logical, check for MvSim.jl package and install if necessary
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @export
bigsimr_setup <- function (pkg_check = TRUE, ...){
  julia <- JuliaCall::julia_setup(installJulia = TRUE, ...)
  JuliaCall::julia_install_package_if_needed("https://github.com/adknudson/MvSim.jl")
  JuliaCall::julia_library("MvSim")
  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(MvSim)),\"!\"=>\"_bang\"))")
  functions <- c(functions, "rand")
  MvSim <- JuliaCall::julia_pkg_import("MvSim", functions)

  MvSim$Pearson  <- JuliaCall::julia_eval("Pearson")
  MvSim$Spearman <- JuliaCall::julia_eval("Spearman")
  MvSim$Kendall  <- JuliaCall::julia_eval("Kendall")

  JuliaCall::autowrap("MvDistribution")

  MvSim
}
