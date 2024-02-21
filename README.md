
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigsimr <a href='https://github.com/SchisslerGroup/r-bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

`bigsimr` is an R package for simulating high-dimensional multivariate
data with a target correlation and arbitrary marginal distributions via
Gaussian copula. It utilizes
[Bigsimr.jl](https://github.com/SchisslerGroup/Bigsimr.jl) for its core
routines. For full documentation and examples, please see the
[Bigsimr.jl docs](https://schisslergroup.github.io/Bigsimr.jl/stable/).

## Features

- **Pearson matching** - employs a matching algorithm (Xiao and
  Zhou 2019) to account for the non-linear transformation in the
  Normal-to-Anything (NORTA) step
- **Spearman and Kendall matching** - Use explicit transformations
  (Lebrun and Dutfoy 2009)
- **Nearest Correlation Matrix** - Calculate the nearest positive
  \[semi\]definite correlation matrix (Qi and Sun 2006)
- **Fast Approximate Correlation Matrix** - Calculate an approximation
  to the nearest positive definite correlation matrix
- **Random Correlation Matrix** - Generate random positive
  \[semi\]definite correlation matrices
- **Fast Multivariate Normal Generation** - Utilize multithreading to
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
bs <- bigsimr_setup()
dist <- distributions_setup()

set.seed(2024-02-20)
```

### Examples

Pearson matching

``` r
(target_corr <- bs$cor_randPD(3))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1770422  0.2197788
#> [2,] -0.1770422  1.0000000 -0.8153085
#> [3,]  0.2197788 -0.8153085  1.0000000

margins <- c(
  dist$Binomial(20, 0.2), 
  dist$Beta(2, 3), 
  dist$LogNormal(3, 1)
)

(adjusted_corr <- bs$pearson_match(target_corr, margins))
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1820291  0.2874494
#> [2,] -0.1820291  1.0000000 -0.9941172
#> [3,]  0.2874494 -0.9941172  1.0000000

x <- bs$rvec(100000, adjusted_corr, margins)

bs$cor(x, bs$Pearson)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1712600  0.2117211
#> [2,] -0.1712600  1.0000000 -0.6679074
#> [3,]  0.2117211 -0.6679074  1.0000000
```

Spearman/Kendall matching

``` r
(spearman_corr <- bs$cor_randPD(3))
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.4768102 0.8865416
#> [2,] 0.4768102 1.0000000 0.5318276
#> [3,] 0.8865416 0.5318276 1.0000000

(adjusted_corr <- bs$cor_convert(spearman_corr, bs$Spearman, bs$Pearson))
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.4941437 0.8954010
#> [2,] 0.4941437 1.0000000 0.5497588
#> [3,] 0.8954010 0.5497588 1.0000000

x <- bs$rvec(100000, adjusted_corr, margins)

bs$cor(x, bs$Spearman)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.4663302 0.8746458
#> [2,] 0.4663302 1.0000000 0.5276638
#> [3,] 0.8746458 0.5276638 1.0000000
```

Nearest correlation matrix

``` r
s <- bs$cor_randPSD(200)
r <- bs$cor_convert(s, bs$Spearman, bs$Pearson)
bs$is_correlation(r)
#> [1] FALSE
```

``` r
p <- bs$cor_nearPD(r)
bs$is_correlation(p)
#> [1] TRUE
```

Fast approximate nearest correlation matrix

``` r
s <- bs$cor_randPSD(2000)
r <- bs$cor_convert(s, bs$Spearman, bs$Pearson)
bs$is_correlation(r)
#> [1] FALSE
```

``` r
p <- bs$cor_fastPD(r)
bs$is_correlation(p)
#> [1] TRUE
```

## Issues

This package is just a wrapper for the Julia package. Please file any
bug reports or feature requests over at the
[Bigsimr.jl](https://github.com/SchisslerGroup/Bigsimr.jl/issues)
package repo.
