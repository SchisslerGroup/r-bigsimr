# Convert uniform to marginal via X := F^-1(U)
.u2m <- function(u, margin) {
  margin$p <- quote(u)
  eval(margin)
}

#' Creates \code{n} observations from a multivariate distribution with the
#' given marginals and correlation matrix.
#'
#' @param n The number random vectors to generate.
#' @param rho The input correlation matrix.
#' @param margins The marginal distributions (Typically R's "quantile functions")
#' @param type The type of correlation matrix that is being passed (Assumed to
#'    be Pearson by default).
#' @param ensure_PSD Ensure that the converted correlation matrix is positive
#'    semi-definite. More import if the input correlation type is Kendall or
#'    Spearman.
#' @param cores The number of cores to run on
#' @return A matrix of random vectors generated from the specified marginals
#'   and parameters.
#' @export
rvec <- function(n, rho, margins, type = c("pearson", "kendall", "spearman"),
                 ensure_PSD=FALSE, cores = 1L){

  type <- match.arg(type)
  rho  <- cor_convert(rho, from = type, to = "pearson")
  d    <- length(margins)

  if (ensure_PSD)
    rho <- cor_nearPSD(rho)

  if (!is.integer(cores)) {
    message("Number of cores implicitly cast as an integer.")
    cores <- as.integer(cores)
  }

  # generate multivariate uniform distribution (via Z -> U)
  U <- .rmvuu(n, rho)

  # Apply the copula algorithm
  d <- nrow(rho)

  if (cores <= 1L) {

    rv <- sapply(1:d, function(i){
      .u2m(U[,i], margins[[i]])
    })


  } else {

    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)

    rv <- foreach::foreach(i = 1:d, .combine = 'cbind') %dopar% {
      .u2m(U[,i], margins[[i]])
    }

    parallel::stopCluster(cl)

  }

  rv
}
