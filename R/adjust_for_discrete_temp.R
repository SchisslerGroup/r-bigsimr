adjustForDiscrete <- function(r_hat, param_list, type = "spearman", numSigmas=5){

  ## store values

  ## tmp_pair <- idx_mat[,2]
  param_dat <- do.call("rbind",
                       lapply(param_list, function(param, sigmas = numSigmas){
    tmp_dat <- as.data.frame(param, stringsAsFactors = FALSE)
    ## compute mean and variance
    tmp_dat$mu <- tmp_dat$size*(1-tmp_dat$prob) / tmp_dat$prob
    tmp_dat$sigma2 <- tmp_dat$size^2*(1-tmp_dat$prob) / tmp_dat$prob
    ## sim_dat <- rnbinom(n = 10e5, size = tmp_dat$size, prob = tmp_dat$prob)
    ## mean(sim_dat) ; var(sim_dat)
    ## compute 5 standard deviations above the mean
    tmp_dat$upper <- ceiling(tmp_dat$mu + sigmas*sqrt(tmp_dat$sigma2))
    return(tmp_dat)
  }))
  ## Find last point to evaluate in the finite sum
  upperBound <- max(param_dat$upper)

  ## add col for the discrete-adjustment factor
  param_dat$adj <- unlist(lapply(param_list, function(param){
    ## only compute if it is unique, otherwise look up
    ## tmp_dist <- param[[1]]
    ## tmp_size <-
    ## approximate the infinite sum
    pmf_x <- do.call(what = paste0("d", param[[1]]),
                     args = c(list(x = 0:upperBound), param[-1]))

    ## compute finite sum of cubed pmf
    (1 - sum(pmf_x^3))^(1/2)
  }))

  ## find all pairs
  idx_mat <- combn(x = length(param_list), m = 2)

  ## adjust r_hat (compute rescaled Spearmans)
  r_hat[lower.tri(r_hat)] <- r_hat[upper.tri(r_hat)] <-
    apply(idx_mat, MARGIN = 2, function(tmp_pair){
      ## lookup adjustment
      adj_factor <- param_dat$adj[tmp_pair[1]] * param_dat$adj[tmp_pair[2]]
      tmp_r <- r_hat[tmp_pair[1], tmp_pair[2]]
      ## adj_factor * tmp_r
      tmp_r / adj_factor
    })

  ## return adjusted matrix
  return(r_hat)
}
