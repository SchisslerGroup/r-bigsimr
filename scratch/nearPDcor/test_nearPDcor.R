devtools::load_all()

rho <- matrix(c(
  0.99, 0.78, 0.59, 0.44,
  0.78, 0.92, 0.28, 0.81,
  0.59, 0.28, 1.12, 0.23,
  0.44, 0.81, 0.23, 0.99
), 4, 4, byrow = TRUE)

rho <- cov2cor(rho)

rho_pd1 <- as.matrix(Matrix::nearPD(rho, corr = TRUE, ensureSymmetry = FALSE)$mat)
rho_pd2 <- nearPDcor(rho)

round(rho, 4)
round(rho_pd1, 4)
