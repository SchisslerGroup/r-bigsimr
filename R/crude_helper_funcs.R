#' Adjust the correlation matrix when there are discrete distributions present
#'
#' @param rho The input correlation matrix
#' @param params The parameters of the marginals.
#' @param nSigmas The number of standard deviations from the mean
#' @export
adjustForDiscrete <- function(rho, params, nSigmas) {
  upper_bound <- lapply(params, function(param) {
    prob <- param[["prob"]]
    size <- param[["size"]]

    mu <- size * (1 - prob) / prob
    sigma2 <- size^2 * (1 - prob) / prob

    ceiling(mu + nSigmas * sqrt(sigma2))
  })
  upper_bound <- max(unlist(upper_bound))

  adj <- lapply(params, function(param) {
    if (param[[1]] %in% discrete_dists) {
      pmf_x <- do.call(paste0("d", param[[1]]),
                       c(list(x = 0:upper_bound), param[-1]))
      return(sqrt(1 - sum(pmf_x^3)))
    } else {
      return(1)
    }
  })
  adj <- unlist(adj)

  idx_mat <- combn(x = length(params), m = 2)

  rho_adj <- apply(idx_mat, 2, function(pair) {
    adj_factor <- adj[pair[1]] * adj[pair[2]]
    rho_tmp <- rho[pair[1], pair[2]]
    min(rho_tmp / adj_factor, 1)
  })
  rho[lower.tri(rho)] <- rho[upper.tri(rho)] <- rho_adj
  rho
}
