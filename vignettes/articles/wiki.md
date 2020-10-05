Welcome to the bigsimr wiki!

## Purpose

## Conventions

### Branches

* **master**
  * Used for stable releases and MINOR releases (see the tags section)
  * Changes to the master branch must be made through a pull request
* **develop**
  * Used for new and potentially unstable features
  * PATCH releases will be incorporated into this branch
  * Any member of the SchisslerGroup or member who is a part of working on bigsimr should have read/write permissions for this branch
* **gh-pages**
  * This is where the online documentation gets uploaded to
  * There is no reason to push to this branch manually

### Tags

* The version number in the description file loosely adheres to [Semantic Versioning](https://semver.org/)
* Until a 1.0 release is agreed upon, there can be backwards compatibility-breaking changes in MINOR and PATCH releases
* Generally we will only release tags for MAJOR and MINOR versions
* PATCH release tags may become more common when there is a 1.0 release

### Style and organization

* This package is built around the conventions described in [R Packages](https://r-pkgs.org/)
* Please use the [tidyverse style guide](https://style.tidyverse.org/) when formatting your code
* The master branch should only contain the code necessary for the user and to pass CRAN checks
  * Extra stuff in the master branch must have a rule in the `.Rbuildignore` file
  * notebooks and scratch code should go in the develop branch or personal branches
* For clarity (and CRAN), external (imported) and non base R functions must be prefixed with their namespace
  * If an external package `PKG` is being used, you must first call `usethis::use_package("PKG")` so that it gets added to the imports in DESCRIPTION
  * Examples of when to use namespaces:
    * `stats::rnorm()`
    * `utils::combn()`
    * `mvnfast::rmvn()`

### Coding Philosophy

* Functions should do one (small) thing and do it well
* Functions should be composable
* Code should be self-documenting
  * Comments should provide _extra_ information that cannot be inferred from the code
  * If the code does not clearly describe what it does, then the code should be re-written so that it does
* Writing good tests now can save a lot of time later (see the Testing section)

#### File and Function Naming

* Principal functions (like `rvec`) should be either self-explanatory or follow a common convention
* Functions that can be grouped together should be prefixed by what they have in common
  * E.g. many of the functions that act on or produce correlations are prefixed with `cor_`
  * This also has utility when it comes to documentation; related functions can be grouped with `starts_with()` (see `./pkgdown/_pkgdown.yml`)
* Each primary function should be placed in its own file
  * Any utility functions that it depends on must be placed in the same file
  * It is recommended to use 'dot' naming for these utility functions (e.g. `.set_omega <- function(...)`) because dot functions are generally not exported for the user
* Closely related functions may be placed in the same file
  * Make sure the file name is indicative of all functions that it contains

### Documentation (via roxygen2)

See [object documentation](https://r-pkgs.org/man.html). `@import` and `@importFrom` statements are required so that the NAMESPACE file gets updated when calling `devtools::document()`.

* Every primary function must have full documentation
  * A description of what the function does (including what it takes in and returns)
  * A description for each argument (`@param`)
  * A description of what the function returns (`@return`)
  * At least one example (`@example`)
  * Any necessary imports (`@import`, `@importFrom`) <- this is for updating the NAMESPACE file
  * An export statement (`@export`)
* Secondary functions can be minimally documented
  * A short description of what the function does
  * A description for each argument (`@param`)
  * Any necessary imports (`@import`, `@importFrom`)
  * An export statement (`@export`)
* Utility (not exported) functions must have import statements (`@import`, `@importFrom`)

### Testing

* Write tests early and check tests often
* [Kevlin Henney](https://www.youtube.com/watch?v=azoucC_fwzw) has a lot to say about good unit testing and communication
* Unit tests must be written in the `./tests/testthat` directory and named clearly
* Tests should check a concept
  * E.g. `cor_randPD()` must generate a correlation matrix. What are the requirements for a matrix to be a correlation matrix?
    * symmetric
    * positive definite
    * ones in the diagonal
    * all values in the range [-1, 1]

## Issues

Please submit issues [here](https://github.com/SchisslerGroup/bigsimr/issues) with the appropriate tags

