#' Computes the theoretical upper and lower bounds of possible correlations
#'   given a set of marginals
#'
#' @param margins The parameters of the marginals.
#' @param type The type of correlation matrix that is being passed.
#' @param cores The number of cores to utilize.
#' @param reps The number of sims used to estimate the bounds.
#' @return A list containing the theoretical upper and lower bounds
#' @export
cor_bounds <- function(margins,
                       type = c("pearson", "kendall", "spearman"),
                       cores = 1,
                       reps = 1e5){

  type <- match.arg(type)
  d <- length(margins)
  index_mat <- utils::combn(x = d, m = 2)

  # Replace the quantile function with the RNG function (e.g. qnorm -> rnorm)
  q2r <- function(x) {
    s <- rlang::as_string(x)
    substr(s, 1, 1) <- "r"    # replace 'q' with 'r'
    rlang::ensym(s)
  }

  # Generate random samples for each margin and sort the vectors
  # The sorted vectors are used for computing the bounds
  if (cores == 1) {
    sim_data <- sapply(1:d, function(i) {
      # Replace the quantile function with the RNG function (e.g. qnorm -> rnorm)
      margins[[i]][[1]] <- q2r(margins[[i]][[1]])
      margins[[i]]$n <- quote(reps)
      eval(rlang::call2("sort", margins[[i]]))
    })

    # Upper bounds
    rho_upper <- cor_fast(sim_data, method = type)

    # Lower bounds
    rho_lower_values <- apply(index_mat, 2, function(index, data, ...){
      cor_fast(data[, index[1]],
               rev(data[, index[2]]),
               method = type)
    }, data = sim_data, method = type)

    rho_lower <- matrix(0, d, d)
    diag(rho_lower) <- 0.5
    rho_lower[lower.tri(rho_lower)] <- rho_lower_values
    rho_lower <- rho_lower + t(rho_lower)

  } else {
    # Set up the parallel computing environment
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)

    sim_data <- foreach::foreach(i = 1:d, .combine = "cbind") %dopar% {
      margins[[i]][[1]] <- q2r(margins[[i]][[1]])
      margins[[i]]$n <- quote(reps)
      eval(rlang::call2("sort", margins[[i]]))
    }

    # Upper bounds
    rho_upper_values <- parallel::parApply(
      cl = cl,
      X = index_mat,
      MARGIN = 2,
      FUN = function(index, data, ...){
        cor_fast(x = data[, index[1]],
                 y = data[, index[2]],
                 method = type)
      }, data = sim_data, method = type)

    rho_upper <- matrix(0, d, d)
    diag(rho_upper) <- 0.5
    rho_upper[lower.tri(rho_upper)] <- rho_upper_values
    rho_upper <- rho_upper + t(rho_upper)

    # Lower bounds
    rho_lower_values <- parallel::parApply(
      cl     = cl,
      X      = index_mat,
      MARGIN = 2,
      FUN    = function(index, data, ...){
        cor_fast(x      = data[, index[1]],
                 y      = rev(data[, index[2]]),
                 method = type)
      }, data = sim_data, method = type)

    rho_lower <- matrix(0, d, d)
    diag(rho_lower) <- 0.5
    rho_lower[lower.tri(rho_lower)] <- rho_lower_values
    rho_lower <- rho_lower + t(rho_lower)

    # Shut down the parallel cluster
    parallel::stopCluster(cl)
  }

  list(upper = rho_upper,
       lower = rho_lower)
}
