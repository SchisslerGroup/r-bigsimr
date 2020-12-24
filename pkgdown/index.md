# bigsimr <a href='https://github.com/SchisslerGroup/bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->
[![Documentation](https://img.shields.io/badge/docs-release-blue.svg)](https://schisslergroup.github.io/bigsimr/reference/index.html)
[![Build
Status](https://travis-ci.com/SchisslerGroup/bigsimr.svg?branch=master)](https://travis-ci.com/SchisslerGroup/bigsimr)
[![Licence](https://img.shields.io/github/license/schisslergroup/bigsimr)](https://choosealicense.com/licenses/gpl-3.0/)
![Release](https://img.shields.io/github/v/tag/schisslergroup/bigsimr?label=release&sort=semver)
<!-- badges: end -->


`bigsimr` is an R package for simulating high-dimensional multivariate data with arbitrary marginal distributions

bigsimr lets you simulate multivariate data given a correlation matrix and a list of distributions. The correlation matrix can be of type Pearson, Spearman, or Kendall.

You can install the release version of the package from GitHub:

```r
devtools::install_github("SchisslerGroup/bigsimr")
```

To get a bug fix or to use a new feature, you can install the development version from GitHub:

```r
devtools::install_github("SchisslerGroup/bigsimr", ref="develop")
```

This package depends on
[reticulate](https://rstudio.github.io/reticulate/) to draw on the speed
of Google’s [jax](https://github.com/google/jax) library. Please see the [bigsimr installation instructions](https://schisslergroup.github.io/bigsimr/articles/install-bigsimr.html) for more details.

If on Windows or a system without Python, then the package will default to alternative methods. The option can be toggled with

```r
options(use_jax = FALSE)
```

Currently Google's jax library does not have ready-to-use binaries for Windows. The recommendation is to use Windows Subsystem for Linux (WSL), otherwise you will need to build `jaxlib` from source ([see here](https://jax.readthedocs.io/en/latest/developer.html#additional-notes-for-building-jaxlib-from-source-on-windows)).

Additionally, Windows does not allow for forked multiprocessing, so there will be no performance enhancement on a multicore Windows machine.
