readCSV <- function(file) {
  r <- read.csv(file)
  r <- as.matrix(r)
  attr(r, 'dimnames') <- NULL
  r[,-1]
}

# Create a negative definite correlation matrix
rho <- matrix(c(
  0.99, 0.78, 0.59, 0.44,
  0.78, 0.92, 0.28, 0.81,
  0.59, 0.28, 1.12, 0.23,
  0.44, 0.81, 0.23, 0.99
), 4, 4, byrow = TRUE)
rho1 <- cov2cor(rho)

rho2 <- readCSV("scratch/rho_ND_1026.csv")
rho3 <- readCSV("scratch/rho_ND_3076.csv")
rho4 <- readCSV("scratch/rho_ND_5127.csv")


system.time( r1_2 <- Matrix::nearPD(rho1, corr = TRUE, ensureSymmetry = FALSE) )
system.time( r2_2 <- Matrix::nearPD(rho2, corr = TRUE, ensureSymmetry = FALSE) )
system.time( r3_2 <- Matrix::nearPD(rho3, corr = TRUE, ensureSymmetry = FALSE) )
]3system.time( r4_2 <- Matrix::nearPD(rho4, corr = TRUE, ensureSymmetry = FALSE) )

