
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigsimr <a href='https://github.com/adknudson/bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

<!-- badges: end -->

Simulate arbitrary multivariate distributions efficiently in R.

## Installation

You can install the released version of bigsimr from
[github](https://github.com/) with:

``` r
# install.packages("devtools")
remotes::install_github("SchisslerGroup/bigsimr")
```

This package relies on
[reticulate](https://rstudio.github.io/reticulate/) to draw on the speed
of Google’s [jax](https://github.com/google/jax) library. The easiest
way to get a working python environment is to use miniconda and create a
new environment. We provide a handy function for making the conda
environment.

``` r
library(bigsimr)
install_bigsimr(method = "conda", envname = "bigsimr")
```

### Development version

To get a bug fix or to use a feature from the development version, you
can install the development version of bigsimr from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("SchisslerGroup/bigsimr", ref="develop")
```

### GPU version

Since jax has the ability to compile code to the CPU and GPU, bigsimr
also has the ability to benefit from GPU acceleration. Follow the [jax
installation guide](https://github.com/google/jax#installation) for
getting a GPU version.

## Overview

  - `rvec()` simulates multivariate data with specified marginal
    distributions and correlation
  - `rcor()` generate a random correlation matrix
  - `convertCor()` convert *from* one correlation *to* another
  - `computeCorBounds()` compute the theoretical lower and upper
    correlations for a target multivariate distribution

## Usage

``` r
library(bigsimr)
# Reticulate needs to be able to find the python binary with `jax` installed
reticulate::use_condaenv("bigsimr-cpu")
```

to generate multivariate data, we need a list of marginals (and their
parameters), and a correlation structure (matrix). The marginal
distributions can be built up using R’s special `alist` function. This
allows one to enter the distributions without evaluating anything (yet).

``` r
margins = alist(
  qnorm(mean = 3.14, sd = 0.1),
  qbeta(shape1 = 1, shape2 = 4),
  qnbinom(size = 10, prob = 0.75)
)
```

The next step is to define a correlation structure for the multivariate
distribution. This correlation matrix can either come from observed
data, or we can set it ourselves, or we can generate a random
correlation matrix via `bigsimr::rcor()`.

``` r
# rho <- rcor(d = 3)

rho <- matrix(0.5, nrow = 3, ncol = 3)
diag(rho) <- 1.0
rho
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.5  0.5
#> [2,]  0.5  1.0  0.5
#> [3,]  0.5  0.5  1.0
```

Finally we can generate a random vector with our specified marginals and
correlation structure. The last argument, `type`, is looking to know
what kind of correlation matrix it is receiving. Right now it can handle
the Pearson product-moment correlation, Spearman’s \(\rho\), or
Kendall’s \(\tau\).

``` r
x <- rvec(100, rho = rho, margins = margins, type = "pearson")
```

``` r
# Sample correlation
cor(x)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.4959944 0.5525602
#> [2,] 0.4959944 1.0000000 0.5496653
#> [3,] 0.5525602 0.5496653 1.0000000
```

``` r
# Estimated upper and lower correlation bounds
computeCorBounds(margins, type = "pearson")
#> $upper
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.9533036 0.9708210
#> [2,] 0.9533036 1.0000000 0.9828064
#> [3,] 0.9708210 0.9828064 1.0000000
#> 
#> $lower
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.9534858 -0.9706092
#> [2,] -0.9534858  1.0000000 -0.8696132
#> [3,] -0.9706092 -0.8696132  1.0000000
```

## Appendix

``` r
all_dists <- alist(
  qbeta( shape1, shape2 ),
  qbinom( size, prob ),
  qcauchy( location, scale ),
  qchisq( df ),
  qexp( rate ),
  qf( df1, df2 ),
  qgamma( shape, rate ),
  qgeom( prob ),
  qhyper( m, n, k ),
  qlogis( location, scale ),
  qlnorm( meanlog, sdlog ),
  qnbinom( size, prob ),
  qnorm( mean, sd ),
  qpois( lambda ),
  qt( df ),
  qunif( min, max ),
  qweibull( shape, scale ),
  qwilcox( m, n ),
  qsignrank( n )
)
```
