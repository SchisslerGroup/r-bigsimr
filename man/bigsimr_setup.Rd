% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigsimr.R
\name{bigsimr_setup}
\alias{bigsimr_setup}
\title{Setup bigsimr}
\usage{
bigsimr_setup(pkg_check = TRUE, ...)
}
\arguments{
\item{pkg_check}{logical, check for Bigsimr.jl package and install if necessary}

\item{...}{Parameters are passed down to JuliaCall::julia_setup}
}
\value{
Return the imported wrapper of Bigsimr.jl Julia package
}
\description{
This function initializes Julia and the Bigsimr.jl package.
The first time will be longer since it includes precompilation.
Additionally, this will install Julia and the required packages
if they are missing.
}
