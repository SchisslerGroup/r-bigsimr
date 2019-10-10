#' Creates \code{n} observations from a multivariate distribution with the
#' given marginals and correlation matrix.
#'
#' @param n The number random vectors to generate.
#' @param R The input correlation matrix.
#' @param params The parameters of the marginals.
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed.
#' @return A matrix of random vectors generated from the specified marginals and parameters.
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
rvec <- function(n, R, params, cores = 1,
                 type = c("pearson", "kendall", "spearman")){
  # Handle different types of dependencies
  if (type == "spearman") {
    R <- convertSpearmanPearson(R)
  }
  if (type == "kendall") {
    R <- convertKendallPearson(R)
  }

  # determine the dimension d
  d <- NROW(R)

  # generate MVN sample
  mvn_sim <- mvnfast::rmvn(n = n,
                           mu = rep(0, d),
                           sigma = R,
                           ncores = cores,
                           isChol = FALSE)

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
  colnames(mv_sim) <- rownames(R)
  return(mv_sim)
}
