
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigsimr <a href='https://github.com/SchisslerGroup/r-bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

`bigsimr` is an R package for simulating high-dimensional multivariate
data with a target correlation and arbitrary marginal distributions via
Gaussian copula. It utilizes
[Bigsimr.jl](https://github.com/adknudson/Bigsimr.jl) for its core
routines. For full documentation and examples, please see the
[Bigsimr.jl docs](https://adknudson.github.io/Bigsimr.jl/stable/).

## Features

-   **Pearson matching** - employs a matching algorithm (Xiao and
    Zhou 2019) to account for the non-linear transformation in the
    Normal-to-Anything (NORTA) step
-   **Spearman and Kendall matching** - Use explicit transformations
    (Lebrun and Dutfoy 2009)
-   **Nearest Correlation Matrix** - Calculate the nearest positive
    [semi]definite correlation matrix (Qi and Sun 2006)
-   **Fast Approximate Correlation Matrix** - Calculate an approximation
    to the nearest positive definite correlation matrix
-   **Random Correlation Matrix** - Generate random positive
    [semi]definite correlation matrices
-   **Fast Multivariate Normal Generation** - Utilize multithreading to
    generate multivariate normal samples in parallel

## Installation

You can install the release version of the package from GitHub:

``` r
remotes::install_github("SchisslerGroup/r-bigsimr")
```

To get a bug fix or to use a new feature, you can install the
development version from GitHub:

``` r
remotes::install_github("SchisslerGroup/r-bigsimr", ref="develop")
```

Note that the first invocation of `bigsimr::bigsimr_setup()` will
install both Julia and the required packages if they are missing. If you
wish to have it use an existing Julia binary, make sure that `julia` is
found in the path. For more information see the `julia_setup()` function
from [JuliaCall](https://github.com/Non-Contradiction/JuliaCall).

## Usage

``` r
library(bigsimr)
Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores()) # activate multithreading
bs <- bigsimr_setup()
dist <- distributions_setup()

set.seed(2020-02-28)
```

### Examples

Pearson matching

``` r
(target_corr <- bs$cor_randPD(3))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.6499761 -0.4650075
#> [2,] -0.6499761  1.0000000  0.4249699
#> [3,] -0.4650075  0.4249699  1.0000000
margins <- c(dist$Binomial(20, 0.2), dist$Beta(2, 3), dist$LogNormal(3, 1))
(adjusted_corr <- bs$pearson_match(target_corr, margins))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.6699159 -0.6484545
#> [2,] -0.6699159  1.0000000  0.5510696
#> [3,] -0.6484545  0.5510696  1.0000000
x <- bs$rvec(100000, adjusted_corr, margins)
bs$cor(x, bs$Pearson)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.6494003 -0.4684822
#> [2,] -0.6494003  1.0000000  0.4301399
#> [3,] -0.4684822  0.4301399  1.0000000
```

Spearman/Kendall matching

``` r
(spearman_corr <- bs$cor_randPD(3))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.7330679 -0.1263919
#> [2,]  0.7330679  1.0000000 -0.1690453
#> [3,] -0.1263919 -0.1690453  1.0000000
(adjusted_corr <- bs$cor_convert(spearman_corr, bs$Spearman, bs$Pearson))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.7489555 -0.1322607
#> [2,]  0.7489555  1.0000000 -0.1767928
#> [3,] -0.1322607 -0.1767928  1.0000000
x <- bs$rvec(100000, adjusted_corr, margins)
bs$cor(x, bs$Spearman)
#>            [,1]      [,2]       [,3]
#> [1,]  1.0000000  0.721592 -0.1261553
#> [2,]  0.7215920  1.000000 -0.1729510
#> [3,] -0.1261553 -0.172951  1.0000000
```

Nearest correlation matrix

``` r
s <- bs$cor_randPSD(200)
r <- bs$cor_convert(s, bs$Spearman, bs$Pearson)
bs$iscorrelation(r)
#> [1] FALSE
```

``` r
p <- bs$cor_nearPD(r)
bs$iscorrelation(p)
#> [1] TRUE
```

Fast approximate nearest correlation matrix

``` r
s = bs$cor_randPSD(2000)
r = bs$cor_convert(s, bs$Spearman, bs$Pearson)
bs$iscorrelation(r)
#> [1] FALSE
```

``` r
p = bs$cor_fastPD(r)
bs$iscorrelation(p)
#> [1] TRUE
```

------------------------------------------------------------------------

<!-- badges: start -->

[![Documentation](https://img.shields.io/badge/docs-release-blue.svg)](https://schisslergroup.github.io/bigsimr/reference/index.html)
[![Build
Status](https://travis-ci.com/SchisslerGroup/bigsimr.svg?branch=master)](https://travis-ci.com/SchisslerGroup/bigsimr)
![Release](https://img.shields.io/github/v/tag/schisslergroup/bigsimr?label=release&sort=semver)
<!-- badges: end -->
