
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bigsimr <a href='https://github.com/SchisslerGroup/bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

### bigsimr is an R package for simulating high-dimensional multivariate data with arbitrary marginal distributions

bigsimr lets you simulate multivariate data given a correlation matrix
and a list of distributions. The correlation matrix can be of type
Pearson, Spearman, or Kendall, and we use a matching algorithm to ensure
that the estimated correlation of the simulated data is the same as the
input correlation.

### See the [website](https://schisslergroup.github.io/bigsimr/) for more information, including [installation instructions](https://schisslergroup.github.io/bigsimr/articles/install-bigsimr.html), [tutorials](https://schisslergroup.github.io/bigsimr/articles/using-rvec.html), and [package documentation](https://schisslergroup.github.io/bigsimr/reference/index.html).

You can install the release version of the package from GitHub:

``` r
devtools::install_github("SchisslerGroup/bigsimr")
```

To get a bug fix or to use a new feature, you can install the
development version from GitHub:

``` r
devtools::install_github("SchisslerGroup/bigsimr", ref="develop")
```

This package depends on
[reticulate](https://rstudio.github.io/reticulate/) to draw on the speed
of Google’s [jax](https://github.com/google/jax) library. Please see the
[bigsimr installation instructions](#) for more details.

-----

<!-- badges: start -->

[![Documentation](https://img.shields.io/badge/docs-release-blue.svg)](https://schisslergroup.github.io/bigsimr/reference/index.html)
[![Build
Status](https://travis-ci.com/SchisslerGroup/bigsimr.svg?branch=master)](https://travis-ci.com/SchisslerGroup/bigsimr)
[![Licence](https://img.shields.io/github/license/schisslergroup/bigsimr)](https://choosealicense.com/licenses/gpl-3.0/)
![Release](https://img.shields.io/github/v/tag/schisslergroup/bigsimr?label=release&sort=semver)
<!-- badges: end -->
