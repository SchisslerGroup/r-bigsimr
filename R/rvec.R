#' Creates \code{n} observations from a multivariate distribution with the
#' given marginals and correlation matrix.
#'
#' @param n The number random vectors to generate.
#' @param rho The input correlation matrix.
#' @param margins The marginal distributions (Typically R's "quantile functions")
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed (Assumed to
#'   be Pearson by default).
#' @return A matrix of random vectors generated from the specified marginals
#'   and parameters.
#' @export
rvec <- function(n, rho, margins, cores = 1,
                 type = c("pearson", "kendall", "spearman")){

  type <- match.arg(type)

  rho <- convertCor(rho, from = type, to = "pearson")

  cores <- as.integer(min(parallel::detectCores(), max(cores, 1)))

  # Convert uniform to marginal via X := F^-1(U)
  u2m <- function(u, margin) {
    margin$p <- quote(u)
    eval(margin)
  }

  # generate multivariate uniform distribution (via Z -> U)
  U <- .rmvuu(n, rho)

  # Apply the copula algorithm
  d <- nrow(rho)
  if (cores == 1) {

    rv <- sapply(1:d, function(i){
      u2m(U[,i], margins[[i]])
    })

  } else {

    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)

    rv <- foreach::foreach(i = 1:d, .combine = 'cbind') %dopar% {
      u2m(U[,i], margins[[i]])
    }

    parallel::stopCluster(cl)

  }

  rv
}
