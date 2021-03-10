
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigsimr <a href='https://github.com/SchisslerGroup/r-bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

`bigsimr` is an R package for simulating high-dimensional multivariate
data with a target correlation and arbitrary marginal distributions via
Gaussian copula. It utilizes
[Bigsimr.jl](https://github.com/SchisslerGroup/Bigsimr.jl) for its core
routines. For full documentation and examples, please see the
[Bigsimr.jl docs](https://SchisslerGroup.github.io/Bigsimr.jl/stable/).

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
#> [1,]  1.0000000 -0.3566613 -0.4686234
#> [2,] -0.3566613  1.0000000  0.6501013
#> [3,] -0.4686234  0.6501013  1.0000000
margins <- c(dist$Binomial(20, 0.2), dist$Beta(2, 3), dist$LogNormal(3, 1))
(adjusted_corr <- bs$pearson_match(target_corr, margins))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.3670106 -0.6537634
#> [2,] -0.3670106  1.0000000  0.8431715
#> [3,] -0.6537634  0.8431715  1.0000000
x <- bs$rvec(100000, adjusted_corr, margins)
bs$cor(x, bs$Pearson)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.3592258 -0.4693409
#> [2,] -0.3592258  1.0000000  0.6483799
#> [3,] -0.4693409  0.6483799  1.0000000
```

Spearman/Kendall matching

``` r
(spearman_corr <- bs$cor_randPD(3))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.5696526 -0.5277236
#> [2,]  0.5696526  1.0000000 -0.1333224
#> [3,] -0.5277236 -0.1333224  1.0000000
(adjusted_corr <- bs$cor_convert(spearman_corr, bs$Spearman, bs$Pearson))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.5877329 -0.5456254
#> [2,]  0.5877329  1.0000000 -0.1395015
#> [3,] -0.5456254 -0.1395015  1.0000000
x <- bs$rvec(100000, adjusted_corr, margins)
bs$cor(x, bs$Spearman)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000  0.5633441 -0.5227178
#> [2,]  0.5633441  1.0000000 -0.1322200
#> [3,] -0.5227178 -0.1322200  1.0000000
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

## Issues

This package is just a wrapper for the Julia package. Please file any
bug reports or feature requests over at the
[Bigsimr.jl](https://github.com/SchisslerGroup/Bigsimr.jl/issues)
package repo.
