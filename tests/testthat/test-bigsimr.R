test_that("Bigsimr.jl loads", {
  skip_on_cran()

  bs   <- bigsimr::bigsimr_setup()
  dist <- bigsimr::distributions_setup()
})
