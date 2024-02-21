test_that("Bigsimr.jl loads", {
  skip_on_cran()

  bs   <- bigsimr::bigsimr_setup()
  dist <- bigsimr::distributions_setup()
})



test_that("Pearson matching is working", {
  skip_on_cran()

  bs   <- bigsimr::bigsimr_setup()
  dist <- bigsimr::distributions_setup()

  JuliaCall::julia_eval('using Random; Random.seed!(1);')
  target_corr <- bs$cor_randPD(3)
  margins <- c(dist$Binomial(20, 0.2), dist$Beta(2, 3), dist$LogNormal(3, 1))
  adjusted_corr <- bs$pearson_match(target_corr, margins)
  x <- bs$rvec(100000, adjusted_corr, margins)

  expect_equal(mean(bs$cor(x, bs$Pearson)), mean(target_corr), tolerance = 0.05)
})

test_that("Spearman matching is working", {
  skip_on_cran()

  bs   <- bigsimr::bigsimr_setup()
  dist <- bigsimr::distributions_setup()

  JuliaCall::julia_eval('using Random; Random.seed!(1);')
  spearman_corr <- bs$cor_randPD(3)
  margins <- c(dist$Binomial(20, 0.2), dist$Beta(2, 3), dist$LogNormal(3, 1))
  adjusted_corr <- bs$cor_convert(spearman_corr, bs$Spearman, bs$Pearson)
  x <- bs$rvec(100000, adjusted_corr, margins)

  expect_equal(mean(bs$cor(x, bs$Spearman)), mean(spearman_corr), tolerance = 0.05)
})
