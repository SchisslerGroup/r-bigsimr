devtools::load_all()

#= D=10, N=1e5, Cores=24
d <- 10
set.seed(12)
rho <- rcor(d = d)
params <- rnbinom_params(d, shape = 100, id_margins = FALSE)
n <- 10
cores <- 1


