#' Creates \code{n} observations from a multivariate negative binomial distribution.
#' 
#' @param n The number random vectors to generate.
#' @param R The input correlation matrix.
#' @param params The parameters of the marginals.
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed.
simInvNB <- function(n, R, params, cores = 1, type = c("pearson", "kendall", "spearman")){
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
  mvn_sim <- mvnfast::rmvn(n = n, mu = rep(0, d), sigma = R, ncores = cores, isChol = FALSE)
  
  # single threaded
  if (cores == 1) {
    # 3. now DIRECTLY TO negative binomial inverse transform to gamma
    nb_sim <- sapply(1:d, function(i){
      qnbinom(p = pnorm(mvn_sim[,i]),
              size = params["alpha",i] , 
              prob = (1 + params["lambda",i])^(-1))
    })
  } else { # mulit-core --- parallelize
    # Start cluster
    cl <- parallel::makeCluster(cores, type = "FORK")
    
    # now DIRECTLY TO transform to negative binomial
    doParallel::registerDoParallel(cl)
    nb_sim <- foreach::foreach(i = 1:d, .combine = 'cbind') %dopar% {
      qnbinom(p = pnorm(mvn_sim[,i]), 
              size = params["alpha",i] , 
              prob = (1 + params["lambda",i])^(-1))
    }
    
    # close cluster
    parallel::stopCluster(cl)
  }
  
  # Return the simulated data set
  colnames(nb_sim) <- rownames(R)
  return(nb_sim)
}