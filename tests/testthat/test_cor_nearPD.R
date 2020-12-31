test_that("cor_nearPD returns a positive definite correlation matrix", {

  # define a negative definite correlation matrix
  p <- matrix(c(
    1.00, 0.82, 0.56, 0.44,
    0.82, 1.00, 0.28, 0.85,
    0.56, 0.28, 1.00, 0.22,
    0.44, 0.85, 0.22, 1.00
  ), 4, 4, byrow=TRUE)

  r <- cor_nearPD(p)
  expect_true(is_valid_correlation(r))

})
