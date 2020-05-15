## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

reticulate::use_condaenv("py37")
devtools::load_all()

## ----setup, eval=FALSE--------------------------------------------------------
#   # Optional step if your python binary is in a miniconda environment
#  reticulate::use_condaenv("py37")
#  library(bigsimr)

## -----------------------------------------------------------------------------
# Generating samples from a normal distribution
# We need three things: 
#   1) The number of samples to generate
#   2) The location of the normal distribution (mean)
#   3) The scale of the normal distribution (standard deviation)

rnorm(n = 16, mean = 5, sd = 1.5)

## -----------------------------------------------------------------------------
margins = list(
  list("norm", mean = 3.14, sd = 0.1),
  list("beta", shape1 = 1, shape2 = 4),
  list("nbinom", size = 10, prob = 0.75)
)

## -----------------------------------------------------------------------------
rho <- matrix(0.5, nrow = 3, ncol = 3)
diag(rho) <- 1.0
rho

## -----------------------------------------------------------------------------
x = rvec(10, rho = rho, params = margins, type = "pearson")

## ---- echo=FALSE--------------------------------------------------------------
warning("warning.warn('No GPU/TPU found, falling back to CPU.')")

## -----------------------------------------------------------------------------
x

## ---- fig.width=7-------------------------------------------------------------
x = rvec(10000, rho = rho, params = margins, type = "pearson")

par(mfrow=c(1,3))
hist(x[,1], breaks = 30, xlab = "", main = "Normal")
hist(x[,2], breaks = 30, xlab = "", main = "Beta")
hist(x[,3], breaks = 30, xlab = "", main = "Negative Binomial")

## -----------------------------------------------------------------------------
cor(x)

## ---- eval=FALSE--------------------------------------------------------------
#  all_dists <- list(
#    list(dist = "beta", shape1, shape2),
#    list(dist = "binom", size, prob),
#    list(dist = "cauchy", location, scale),
#    list(dist = "chisq", df),
#    list(dist = "exp", rate),
#    list(dist = "f", df1, df2),
#    list(dist = "gamma", shape, rate),
#    list(dist = "geom", prob),
#    list(dist = "hyper", m, n, k),
#    list(dist = "logis", location, scale),
#    list(dist = "lnorm", meanlog, sdlog),
#    list(dist = "nbinom", size, prob),
#    list(dist = "norm", mean, sd),
#    list(dist = "pois", lambda),
#    list(dist = "t", df),
#    list(dist = "unif", min, max),
#    list(dist = "weibull", shape, scale),
#    list(dist = "wilcox", m, n),
#    list(dist = "signrank", n)
#  )

