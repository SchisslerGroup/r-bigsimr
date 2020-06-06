devtools::load_all()

d <- 5
rho <- rcor(d)
margins <- list(
  list("nbinom", size = 5, prob = 0.3),
  list("exp", rate = 4),
  list("binom", size = 5, prob = 0.7),
  list("norm", mean = 10, sd = 3),
  list("pois", lambda = 10)
)

margins2 <- list(
  list("nbinom", 5, 0.3),
  list("exp", 4),
  list("binom", 5, 0.7),
  list("norm", 10, 3),
  list("pois", 10)
)


cores = parallel::detectCores() - 1
reps = 1e3
type = "spearman"

rho_bounds <- computeCorBounds(margins,
                               cores = cores,
                               type = type,
                               reps = reps)

pm1_rho <- matrix(sample(c(-1, 1), d*d, TRUE), d, d)
diag(pm1_rho) <- 1

constrainRho(rho, rho_bounds)
constrainRho(pm1_rho, rho_bounds)

all_corInBounds(rho, margins, cores, type, rho_bounds)
all_corInBounds(pm1_rho, margins, cores, type, rho_bounds)

which_corInBounds(pm1_rho, rho_bounds)

rvec(10, rho, margins, type = "spearman")
rvec(10, rho, margins, type = "spearman", adjustForDiscrete = TRUE)
