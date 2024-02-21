#' Setup bigsimr
#'
#' This function initializes Julia and the Bigsimr.jl package.
#' The first time will be longer since it includes precompilation.
#' Additionally, this will install Julia and the required packages
#' if they are missing.
#'
#' @param pkg_check logical, check for Bigsimr.jl package and install if necessary
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @return Return the imported wrapper of Bigsimr.jl Julia package
#' @export
bigsimr_setup <- function (pkg_check = TRUE, ...) {
  julia <- JuliaCall::julia_setup(installJulia = TRUE, ...)

  if (pkg_check) {
    JuliaCall::julia_install_package_if_needed("Bigsimr")
  }

  JuliaCall::julia_library("Bigsimr")

  functions <- JuliaCall::julia_eval(
    "filter(isascii, string.(propertynames(Bigsimr)))"
  )

  rm_funcs <- c("Bigsimr",
                "NearestCorrelationMatrix",
                "PearsonCorrelationMatch",
                "CorType",
                "Pearson",
                "Spearman",
                "Kendall")

  functions <- functions[!(functions %in% rm_funcs)]

  bs <- JuliaCall::julia_pkg_import("Bigsimr", functions)
  bs$Pearson  <- JuliaCall::julia_eval("Pearson")
  bs$Spearman <- JuliaCall::julia_eval("Spearman")
  bs$Kendall  <- JuliaCall::julia_eval("Kendall")

  bs
}
