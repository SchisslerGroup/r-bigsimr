source("scratch/nearPDcor/prototype_nearPDcor_utils.R")

readCSV <- function(file) {
  r <- read.csv(file)
  r <- as.matrix(r)
  attr(r, 'dimnames') <- NULL
  r[,-1]
}


nearPDcor <- function(G
  ,tau=1e-5
  ,iter_outer=200
  ,iter_inner=20
  ,N=200
  ,epsilon=1e-2  # pre-cg error tolerance
  ,delta=1e-6    # error tolerance
  ,sigma=1e-4    # Newton method line search tolerance
){
  n <- nrow(G)

  # Make G symmetric and diagonals one
  G <- (G + t(G)) / 2

  b <- rep(1, n)
  if (tau > 0) {
    b = b - tau
    G = G - diag(tau, n, n)
  }
  b0 <- b

  y      <- rep(0, n)
  X      <- G + diag(y)
  eig    <- eigen(X, symmetric = TRUE)
  lambda <- eig$values
  P      <- eig$vectors

  ret   <- R_gradient(y, lambda, P, b0)
  f0    <- ret$f
  Fy    <- ret$Fy
  f     <- f0
  b     <- b0 - Fy

  k       <- 0L
  Omega12 <- set_omega_mat(lambda)
  x0      <- y

  X        <- R_PCA(X, lambda, P)
  val_G    <- 0.5 * norm(G, '2')^2
  val_dual <- val_G - f0
  val_obj  <- 0.5 * norm(X - G, '2')^2
  gap      <- (val_obj - val_dual) / (1 + abs(val_dual) + abs(val_obj))

  normb    <- norm(b, '2')
  normb0   <- norm(b0, '2') + 1
  delta_nb <- normb / normb0

  while ( (gap > delta) && (delta_nb > delta) && (k < iter_outer) ) {

    cv <- R_precond_matrix(Omega12, P)
    d  <- R_pre_cg(b, cv, Omega12, P, epsilon, N)

    slope <- sum((Fy - b0) * d)

    y      <- x0 + d
    X      <- G + diag(y)
    eig    <- eigen(X, symmetric = TRUE)
    lambda <- eig$values
    P      <- eig$vectors
    ret    <- R_gradient(y, lambda, P, b0)
    f      <- ret$f
    Fy     <- ret$Fy

    k_inner <- 0L
    while ( (k_inner <= iter_inner) && (f > f0 + sigma*0.5^k_inner*slope + 1e-6) ) {
      k_inner <- k_inner + 1
      y       <- x0 + d*0.5^k_inner
      X       <- G + diag(y)
      eig     <- eigen(X, symmetric = TRUE)
      lambda  <- eig$values
      P       <- eig$vectors
      ret     <- R_gradient(y, lambda, P, b0)
      f       <- ret$f
      Fy      <- ret$Fy
    }

    x0 <- y
    f0 <- f

    X        <- R_PCA(X, lambda, P)
    val_dual <- val_G - f0
    val_obj  <- 0.5 * norm(X - G, '2')^2
    gap      <- (val_obj - val_dual) / (1 + abs(val_dual) + abs(val_obj))

    b        <- b0 - Fy
    normb    <- norm(b, '2')
    delta_nb <- normb / normb0

    Omega12 <- set_omega_mat(lambda)

    k <- k + 1
  }

  X + diag(tau, n, n)
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


system.time( r1_2 <- nearPDcor(rho1) )
system.time( r2_2 <- nearPDcor(rho2) )
system.time( r3_2 <- nearPDcor(rho3) )
system.time( r4_2 <- nearPDcor(rho4) )
