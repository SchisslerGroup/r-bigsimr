library(fastpipe)

julia <- JuliaCall::julia_setup(JULIA_HOME = "/home/alex/julia-1.5.3/bin/")
JuliaCall::julia_install_package_if_needed("https://github.com/adknudson/MvSim.jl")
JuliaCall::julia_library("MvSim")
funs <- list(
  "cor_nearPD",
  "cor_fastPD",
  "cor_randPD",
  "cor_randPSD",
  "cor_bounds",
  "cor_constrain",
  "cor_convert",
  "cov2cor"
)
MvSim <- JuliaCall::julia_pkg_import("MvSim", func_list = funs)
JuliaCall::autowrap("MvSim.Correlation")
JuliaCall::autowrap("MvSim.Pearson")
JuliaCall::autowrap("MvSim.Spearman")
JuliaCall::autowrap("MvSim.Kendall")



MvSim$cor_randPD(10L)
