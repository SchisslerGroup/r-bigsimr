#' Setup bigsimr
#'
#' This function initializes Julia and the Bigsimr.jl package.
#' The first time will be long since it includes precompilation.
#' Additionally, this will install Julia and the required packages
#' if they are missing.
#'
#' @param pkg_check logical, check for Bigsimr.jl package and install if necessary
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @return Return the imported wrapper of Bigsimr.jl Julia package
#' @examples
#' ## bigsimr_setup() is time-consuming and requires Julia+Bigsimr.jl
#' if (Sys.which('julia')!=''){
#'   library(bigsimr)
#'   bs   <- bigsimr::bigsimr_setup()
#'   dist <- bigsimr::distributions_setup()
#'
#'   JuliaCall::julia_eval('using Random; Random.seed!(1);')
#'   # Generate random target correlation matrix
#'   target_corr <- bs$cor_randPD(3)
#'   # Set the margins of variables
#'   margins <- c(dist$Binomial(20, 0.2), dist$Beta(2, 3), dist$LogNormal(3, 1))
#'   # Adjust target correlation matrix using Pearson matching
#'   adjusted_corr <- bs$pearson_match(target_corr, margins)
#'   # Generate random vectors
#'   x <- bs$rvec(10000, adjusted_corr, margins)
#' }
#' @export
bigsimr_setup <- function (pkg_check = TRUE, ...){
  julia <- JuliaCall::julia_setup(installJulia = TRUE, ...)

  if (pkg_check) JuliaCall::julia_install_package_if_needed("Bigsimr")
  JuliaCall::julia_library("Bigsimr")

  JuliaCall::julia_install_package_if_needed("Distributions")

  functions <- JuliaCall::julia_eval(
    "filter(isascii, replace.(string.(propertynames(Bigsimr)),\"!\"=>\"_bang\"))"
  )
  # TODO: Deprecate importing of specific functions (i.e. "iscorrelation") and
  #       Export all necessary functions in Bigsimr.jl
  functions <- c(functions, "iscorrelation")
  bs <- JuliaCall::julia_pkg_import("Bigsimr", functions)

  bs$Pearson  <- JuliaCall::julia_eval("Pearson")
  bs$Spearman <- JuliaCall::julia_eval("Spearman")
  bs$Kendall  <- JuliaCall::julia_eval("Kendall")

  bs
}
