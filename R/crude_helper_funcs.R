#' Generate a random Pearson correlation matrix
#'
#' @param d A positive integer.
#' @param constant_rho A boolean value.
#' @return A \code{dxd} symmetric matrix with ones on the diagonal.
#' @examples
#' rcor(4)
#' rcor(10, TRUE)
#' @export
rcor <- function(d, constant_rho = FALSE) {
  if (constant_rho) {
    tmp_rho <- runif(n = 1, min = 0.25, max = 0.9)
    sigma <- matrix(data = tmp_rho, nrow = d, ncol = d)
    diag(sigma) <- 1
  } else {
    tmp <- matrix(rnorm(d^2), d, d)
    mcov <- tcrossprod(tmp, tmp)
    sigma <- cov2cor(mcov)
  }
  rownames(sigma) <- colnames(sigma) <- paste0("Var", 1:d)
  return(sigma)
}


#' Generate negative binomial marginal parameters
#'
#' @param d A positive integer.
#' @param id_margins A boolean value specifying if the margins should be identical.
#' @return A matrix with columns containing (in order) lambda, alpha, mean, and variance.
#' @examples
#' rnbinom_params(4)
#' rnbinom_params(10)
#' rnbinom_params(7, TRUE)
rnbinom_params <- function(d, id_margins = FALSE) {
  if (!id_margins) {
    probs <- runif(d, min = 0.3, max = 0.7)
    sizes <- sample(3:8, d, TRUE)
  } else {
    probs <- rep(runif(1, 0.3, 0.7), d)
    sizes <- rep(sample(3:8, 1, TRUE), d)
  }
  purrr::transpose(list(rep("nbinom", d), size = sizes, prob = probs))
}


#' @param rho The input correlation matrix
#' @param params The parameters of the marginals.
#' @param sigma The number of standard deviations from the mean
adjustForDiscrete <- function(rho, params, sigma = 5) {
  upper_bound <- lapply(params, function(param) {
    prob <- param[["prob"]]
    size <- param[["size"]]

    mu <- size * (1 - prob) / prob
    sigma2 <- size^2 * (1 - prob) / prob

    ceiling(mu + sigma * sqrt(sigma2))
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
    rho_tmp / adj_factor
  })
  rho[lower.tri(rho)] <- rho[upper.tri(rho)] <- rho_adj
  rho
}
