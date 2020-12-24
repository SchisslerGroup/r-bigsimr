test_that("cor_nearPSD returns a positive semidefinite correlation matrix", {

  # define a negative definite correlation matrix
  p <- matrix(c(
    1.00, 0.82, 0.56, 0.44,
    0.82, 1.00, 0.28, 0.85,
    0.56, 0.28, 1.00, 0.22,
    0.44, 0.85, 0.22, 1.00
  ), 4, 4, byrow=TRUE)

  r <- cor_nearPSD(p)
  e <- eigen(r)

  # Correlation matrices must:
  # 1. be symmetric
  # 2. have ones on diagonal
  # 3. be positive semi-definite
  # 4. have all values in the domain [-1, 1]
  expect_true(Matrix::isSymmetric(r))
  expect_true(all(diag(r) == 1))
  expect_true(all(e$values >= 0))
  expect_true(all((-1 <= r) & (r <= 1)))

})
