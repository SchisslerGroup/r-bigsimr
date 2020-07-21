test_that("cor_randPD and cor_randPSD generate random correlation matrices", {
  # Correlation matrices must be:
  # 1. symmetric
  # 2. have ones on diagonal
  # 3. be positive semidefinite
  # 4. have all values in the domain [-1, 1]

  d <- 10L

  rho_PD <- cor_randPD(d)
  e <- eigen(rho_PD)
  expect_true(Matrix::isSymmetric(rho_PD))
  expect_true(all(diag(rho_PD) == 1))
  expect_true(all(e$values > 0))
  expect_true(all((-1 <= rho_PD) & (rho_PD <= 1)))

  rho_PSD <- cor_randPSD(d)
  expect_true(Matrix::isSymmetric(rho_PSD))
  expect_true(all(diag(rho_PSD) == 1))
  expect_true(all(e$values >= 0))
  expect_true(all((-1 <= rho_PSD) & (rho_PSD <= 1)))
})
