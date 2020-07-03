source("test/prototype_nearPDcor_utils.R")

# Create a negative definite correlation matrix
rho <- matrix(c(
  0.99, 0.78, 0.59, 0.44,
  0.78, 0.92, 0.28, 0.81,
  0.59, 0.28, 1.12, 0.23,
  0.44, 0.81, 0.23, 0.99
), 4, 4, byrow = TRUE)
rho <- cov2cor(rho)

G <- rho
n <- nrow(G)

# Make G symmetric and diagonals one
G <- (G + t(G)) / 2

y      <- rep(0, n)
X      <- G + diag(y)
eig    <- eigen(X, symmetric = TRUE)
lambda <- eig$values
P      <- eig$vectors

b     <- rep(1, n)
b0    <- b
ret   <- R_gradient(y, lambda, P, b0)
f0    <- ret$f
Fy    <- ret$Fy
f     <- f0
b     <- b0 - Fy
normb <- norm(b, '2')

k          <- 0L
iter_outer <- 200L
iter_inner <- 20L
maxit      <- 200L
tol        <- 1e-2
err_tol    <- 1e-6
sigma1     <- 1e-4

Omega12 <- set_omega_mat(lambda)
x0      <- y

while ( (normb > err_tol) && (k < iter_outer) ) {

  cv <- R_precond_matrix(Omega12, P)
  d  <- R_pre_cg(b, cv, Omega12, P, tol, maxit)

  slope <- sum((Fy - b0) * d)

  y      <- x0 + d
  X      <- G + diag(y)
  eig    <- eigen(X, symmetric = TRUE)
  lambda <- eig$values
  P      <- eig$vectors
  ret    <- R_gradient(y, lambda, P, b0)
  f      <- ret$f
  Fy     <- ret$Fy

  k_inner <- 1L
  while ( (k_inner <= iter_inner) && (f > f0 + sigma1*0.5^k_inner*slope + 1e-5) ) {
    y      <- x0 + d*0.5^k_inner
    X      <- G + diag(y)
    eig    <- eigen(X, symmetric = TRUE)
    lambda <- eig$values
    P      <- eig$vectors
    ret    <- R_gradient(y, lambda, P, b0)
    f      <- ret$f
    Fy     <- ret$Fy

    k_inner <- k_inner + 1
  }

  x0    <- y
  f0    <- f
  b     <- b0 - Fy
  normb <- norm(b, '2')

  Omega12 <- set_omega_mat(lambda)

  k <- k + 1
}

r <- sum(lambda > 0) # number of positive eigenvalues
s <- n - r

rho1 <- cov2cor(R_PCA(X, lambda, P))
rho2 <- as.matrix(Matrix::nearPD(rho, corr = TRUE, ensureSymmetry = FALSE)$mat)

all.equal(rho, rho1)
all.equal(rho, rho2, check.attributes = FALSE)
all.equal(rho1, rho2)

eigen(rho1)
