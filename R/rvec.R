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
#' @examples
#' d <- 5
#' rho <- rcor(d)
#' margins <- alist(
#'   "nbinom", size = 5, prob = 0.3),
#'   list("exp", rate = 4),
#'   list("binom", size = 5, prob = 0.7),
#'   list("norm", mean = 10, sd = 3),
#'   list("pois", lambda = 10)
#' )
#'
#' margins2 <- list(
#'   list("nbinom", 5, 0.3),
#'   list("exp", 4),
#'   list("binom", 5, 0.7),
#'   list("norm", 10, 3),
#'   list("pois", 10)
#' )
#'
#' rvec(10, rho, margins, cores = 1, type = "pearson")
#' rvec(10, rho, margins2, cores = 1, type = "pearson")
#' @export
rvec <- function(n, rho, margins, cores = 1,
                 type = c("pearson", "kendall", "spearman")){

  type <- match.arg(type)
  rho <- convertCor(rho, from = type, to = "pearson")

  # generate multivariate uniform distribution (via Z -> U)
  U <- .rmvuu(n, rho)

  # Apply the NORTA algorithm
  d <- NROW(rho)
  if (cores == 1) {
    rv <- sapply(1:d, function(i){.u2m(U[,i], margins[[i]])})
  } else {
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)
    rv <- foreach::foreach(i = 1:d, .combine = 'cbind') %dopar% {
      .u2m(U[,i], margins[[i]])
    }
    parallel::stopCluster(cl)
  }

  colnames(rv) <- rownames(rho)
  rv
}
