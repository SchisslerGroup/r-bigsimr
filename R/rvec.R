#' Creates \code{n} observations from a multivariate distribution with the
#' given marginals and correlation matrix.
#'
#' @param n The number random vectors to generate.
#' @param rho The input correlation matrix.
#' @param params The parameters of the marginals.
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed.
#' @param adjustForDiscrete A boolean for whether to adjust the correlation
#'   matrix when in the presence of discrete distributions
#' @param nSigmas The number of standard deviations above the mean to set the
#'   upper bound when adjusting for discrete distributions
#' @return A matrix of random vectors generated from the specified marginals
#'   and parameters.
#' @examples
#' d <- 5
#' rho <- rcor(d)
#' margins <- list(
#'   list("nbinom", size = 5, prob = 0.3),
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
rvec <- function(n,
                 rho,
                 params,
                 cores = 1,
                 type = c("pearson", "kendall", "spearman"),
                 adjustForDiscrete = TRUE,
                 nSigmas = 10){
  # Handle different types of dependencies
  if (type == "spearman") {
    my_dists <- unlist(lapply(params, '[[', 1))
    if (adjustForDiscrete && any(my_dists %in% discrete_dists)) {
      rho <- adjustForDiscrete(rho, params, nSigmas)
    }
    rho <- convertCor(rho, from = "spearman", to = "pearson")
  }

  if (type == "kendall") {
    rho <- convertCor(rho, from = "kendall", to = "pearson")
  }

  # determine the dimension d
  d <- NROW(rho)

  # generate MVN sample
  # mvn_sim <- mvnfast::rmvn(n = n,
  #                          mu = rep(0, d),
  #                          sigma = rho,
  #                          ncores = cores,
  #                          isChol = FALSE)
  mvn_sim <- .rmvn_jax(NULL, rep(0, d), rho, as.integer(n))

  if (cores == 1) {
    mv_sim <- sapply(1:d, function(i){
      normal2marginal(mvn_sim[,i], params[[i]])
    })
  } else {
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)
    mv_sim <- foreach::foreach(i = 1:d, .combine = 'cbind') %dopar% {
      normal2marginal(mvn_sim[,i], params[[i]])
    }
    parallel::stopCluster(cl)
  }
  colnames(mv_sim) <- rownames(rho)
  return(mv_sim)
}
