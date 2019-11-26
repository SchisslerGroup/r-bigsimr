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

  rho[lower.tri(rho)] <- rho_adj
  rho[upper.tri(rho)] <- t(rho)[upper.tri(rho)]
  rho
}


#' Computes the theoretical upper and lower bounds of possible correlations
#'   given a set of marginals
#'
#' @param params The parameters of the marginals.
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed.
#' @return A list containing the theoretical upper and lower bounds
#' @export
computeCorBounds <- function(params,
                             cores = 1,
                             type = c("pearson", "kendall", "spearman"),
                             reps = 1e5){

  d <- length(params)

  # Generate random samples for each margin and sort the vectors
  # The sorted vectors are used for computing the bounds
  if (cores == 1) {
    sim_data <- sapply(1:d, function(i) {
      sort(do.call(what = paste0("r", unlist(params[[i]][1])),
                   args = c(list(n = reps), params[[i]][-1])))
    })
  } else {
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)
    sim_data <- foreach::foreach(i = 1:d, .combine = "cbind") %dopar% {
        sort(do.call(what = paste0("r", unlist(params[[i]][1])),
                     args = c(list(n = reps), params[[i]][-1])))
      }
  }

  unname(sim_data)
  index_mat <- combn(x = d, m = 2)


  if (cores == 1) {
    # Upper bounds
    rho_upper <- cor(sim_data, method = type)

    # Lower bounds
    rho_lower_values <- apply(index_mat, 2, function(index, data, ...){
      cor(data[, index[1]],
          rev(data[, index[2]]),
          method = type)
    }, data = sim_data, method = type)

    rho_lower <- diag(x = 1, nrow = d, ncol = d)
    rho_lower[lower.tri(rho_lower)] <- rho_lower_values
    rho_lower <- t(rho_lower)
    rho_lower[lower.tri(rho_lower)] <- rho_lower_values

  } else {
    # Upper bounds
    rho_upper_values <- parallel::parApply(
      cl = cl,
      X = index_mat,
      MARGIN = 2,
      FUN = function(index, data, ...){
        cor(x = data[, index[1]],
            y = data[, index[2]],
            method = type)
    }, data = sim_data, method = type)

    # Ensure that the result is a valid correlation matrix
    rho_upper <- diag(x = 1, nrow = d, ncol = d)
    rho_upper[lower.tri(rho_upper)] <- rho_upper_values
    rho_upper <- t(rho_upper)
    rho_upper[lower.tri(rho_upper)] <- rho_upper_values

    # Lower bounds
    rho_lower_values <- parallel::parApply(
      cl = cl,
      X = index_mat,
      MARGIN = 2,
      FUN = function(index, data, ...){
        cor(x = data[, index[1]],
            y = rev(data[, index[2]]),
            method=type)
      }, data = sim_data, method = type)
    parallel::stopCluster(cl)

    # Ensure that the result is a valid correlation matrix
    rho_lower <- diag(x = 1, nrow = d, ncol = d)
    rho_lower[lower.tri(rho_lower)] <- rho_lower_values
    rho_lower <- t(rho_lower)
    rho_lower[lower.tri(rho_lower)] <- rho_lower_values
  }

  list(upper = rho_upper,
       lower = rho_lower)
}


#' Constrains a correlation matrix to the theoretical upper and lower bounds
#'
#' @param rho The input correlation matrix.
#' @param rho_bounds A list containing the theoretical upper and lower bounds
#' @return The constrained correlation matrix
constrainRho <- function(rho, rho_bounds){
  rho_tmp <- pmin(rho_bounds[["upper"]], rho)
  pmax(rho_bounds[["lower"]], rho_tmp)
}


#' Returns a logical matrix of which correlations are in the feasible region
#'
#' @param rho The input correlation matrix.
#' @param rho_bounds A list containing the theoretical upper and lower bounds
#' @param negate Should the logical values be negated in order to identify
#'   values that are outside the feasible region.
#' @export
which.corInBounds <- function(rho, rho_bounds, negate = FALSE){

  tooSmall <- rho < rho_bounds[["lower"]]
  tooLarge <- rho > rho_bounds[["upper"]]

  outOfBounds <- tooSmall | tooLarge
  inBounds <- !outOfBounds

  # Return either the in bounds or out of bounds matrix
  if (negate) {
    outOfBounds
  } else {
    inBounds
  }
}


#' Returns TRUE if all correlations pairs are within the theoretical bounds
#'
#' @param rho The input correlation matrix.
#' @param params The parameters of the marginals.
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed.
#' @param rho_bounds Pre-computed upper and lower correlation bounds
#' ... Other arguments passed to `computeCoreBounds()`
#' @return Logical. TRUE if all correlations pairs are within the theoretical
#'   bounds, and false otherwise.
all.corInBounds <- function(rho,
                          params,
                          cores = 1,
                          type = c("pearson", "spearman", "kendall"),
                          rho_bounds = NULL, ...){

  if ( is.null(rho_bounds) ){
    rho_bounds <- computeCorBounds(params = params,
                                   cores = cores,
                                   type = type,
                                   ...)
  }

  outOfBounds <- which.corInBounds(rho, rho_bounds, negate = TRUE)

  if (any(outOfBounds)) {
    warning(paste0("At least one of the specified correlations are either ",
                   "too large or too small. Please use 'which.corInBounds() ",
                   "to see which correlations are valid."))
    return(FALSE)
  }
  TRUE
}


#' Estimate the Spearman rank correlation
#'
#' For a Nx2 matrix when N=1e5, fast is about 2x faster, and when N=1e6,
#'   fast is about 3x faster.
estimateSpearmanRho <- function(x, fast = TRUE){
  if (fast) {
    coop::pcor(apply(x, 2, fastrank::fastrank_num_avg))
  } else {
    cor(apply(x, 2, rank))
  }
}
