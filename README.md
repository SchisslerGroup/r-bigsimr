
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
remotes::install_github("adknudson/bigsimr")
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
devtools::install_github("adknudson/bigsimr", ref="develop")
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
reticulate::use_condaenv("bigsimr")
```

To generate multivariate data, we need a list of marginals (and their
parameters), and a correlation structure (matrix). The marginal
distributions can be built up as a list of lists, where each sublist
contains the information for the target distribution.

Note that in each sublist, the first item is an unnamed character string
with the R name of the distribution *without a letter prefix*. E.g.
instead of `rnorm`, we pass in just `"norm"`. The second thing to note
is that the remaining items are *named* arguments that go along with the
distribution. A full list of built-in distributions is found in the
appendix.

``` r
margins = list(
  list("norm", mean = 3.14, sd = 0.1),
  list("beta", shape1 = 1, shape2 = 4),
  list("nbinom", size = 10, prob = 0.75)
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
x <- rvec(100, rho = rho, params = margins, type = "pearson")
```

``` r
# Sample correlation
cor(x)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.3919804 0.4317084
#> [2,] 0.3919804 1.0000000 0.4678965
#> [3,] 0.4317084 0.4678965 1.0000000
```

``` r
# Estimated upper and lower correlation bounds
computeCorBounds(margins)
#> Warning in if (method == "pearson") {: the condition has length > 1 and only the
#> first element will be used

#> Warning in if (method == "pearson") {: the condition has length > 1 and only the
#> first element will be used

#> Warning in if (method == "pearson") {: the condition has length > 1 and only the
#> first element will be used

#> Warning in if (method == "pearson") {: the condition has length > 1 and only the
#> first element will be used
#> $upper
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.9539526 0.9712235
#> [2,] 0.9539526 1.0000000 0.9828965
#> [3,] 0.9712235 0.9828965 1.0000000
#> 
#> $lower
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.9528902 -0.9704918
#> [2,] -0.9528902  1.0000000 -0.8696719
#> [3,] -0.9704918 -0.8696719  1.0000000
```
