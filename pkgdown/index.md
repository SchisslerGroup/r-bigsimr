# bigsimr <a href='https://github.com/SchisslerGroup/bigsimr'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->
[![Documentation](https://img.shields.io/badge/docs-release-blue.svg)](https://schisslergroup.github.io/bigsimr/reference/index.html)
[![Build
Status](https://travis-ci.com/SchisslerGroup/bigsimr.svg?branch=master)](https://travis-ci.com/SchisslerGroup/bigsimr)
[![Licence](https://img.shields.io/github/license/schisslergroup/bigsimr)](https://choosealicense.com/licenses/gpl-3.0/)
![Release](https://img.shields.io/github/v/tag/schisslergroup/bigsimr?label=release&sort=semver)
<!-- badges: end -->


`bigsimr` is an R package for simulating high-dimensional multivariate data with arbitrary marginal distributions

bigsimr lets you simulate multivariate data given a correlation matrix and a list of distributions. The correlation matrix can be of type Pearson, Spearman, or Kendall, and we use a matching algorithm to ensure that the estimated correlation of the simulated data is the same as the input correlation.

