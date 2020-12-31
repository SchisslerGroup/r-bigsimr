# Convert uniform to marginal via X := F^-1(U)
.u2m <- function(u, margin) {
  margin$p <- quote(u)
  eval(margin)
}

#' Simulate correlated multivariate data
#'
#' Creates \code{n} observations from a multivariate distribution with the
#'   given marginals and correlation matrix. The correlation matrix is always
#'   assumed to be a Pearson correlation matrix
#'
#' @param n The number random vectors to generate.
#' @param rho The input (Pearson) correlation matrix.
#' @param margins The marginal distributions (Typically R's "quantile functions")
#' @param cores The number of cores to run on.
#' @param ensure_PSD Ensure that the converted correlation matrix is positive
#'    semi-definite.
#' @return A matrix of random vectors generated from the specified marginals
#'   and parameters.
#' @export
rvec <- function(n, rho, margins, cores = 1L, ensure_PSD = FALSE) {

  d <- length(margins)
  rho_names <- colnames(rho) # save dimnames for rv colnames

  if (ensure_PSD)
    rho <- cor_nearPD(rho, 1e-10)

  if (!is.integer(cores)) {
    message("Number of cores implicitly cast as an integer.")
    cores <- as.integer(cores)
  }

  # generate multivariate uniform distribution (via Z -> U)
  if (getOption("use_jax", FALSE)) {
    U <- .rmvuu(n, rho)
  } else {
    U <- stats::pnorm(mvnfast::rmvn(n, rep(0, d), rho, cores))
  }

  # Apply the copula algorithm
  d <- nrow(rho)

  if (cores <= 1L) {

    rv <- vapply(1:d, function(i){
      .u2m(U[,i], margins[[i]])
    }, numeric(n))


  } else {

    cl <- parallel::makeCluster(spec = cores, type = "FORK")

    rv <- parallel::parSapply(cl = cl, 1:d, function(i) {
      .u2m(U[,i], margins[[i]])
    }, simplify = TRUE)

    parallel::stopCluster(cl)

  }

  colnames(rv) <- rho_names
  rv

}
