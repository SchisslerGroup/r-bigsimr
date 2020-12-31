#' Computes the theoretical upper and lower bounds of possible correlations
#'   given a set of marginals
#'
#' @param margins The parameters of the marginals.
#' @param method The type of correlation matrix that is being passed.
#' @param cores NOT YET IMPLEMENTED
#' @param reps The number of sims used to estimate the bounds.
#' @return A list containing the theoretical upper and lower bounds
#' @export
cor_bounds <- function(margins,
                       method = c("pearson", "kendall", "spearman"),
                       cores = 1,
                       reps = 1e5){

  method <- match.arg(method)
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
  sim_data <- sapply(1:d, function(i) {
    # Replace the quantile function with the RNG function (e.g. qnorm -> rnorm)
    margins[[i]][[1]] <- q2r(margins[[i]][[1]])
    # Add the number of reps as an argument
    margins[[i]]$n <- quote(reps)
    # the below statement equates to: sort(rdist(n, params...))
    eval(rlang::call2("sort", margins[[i]]))
  })

  # Upper bounds
  rho_upper <- cor_fast(sim_data, method = method)

  # Lower bounds
  rho_lower_values <- apply(index_mat, 2, function(index, data, ...){
    cor_fast(data[, index[1]],
             rev(data[, index[2]]),
             method = method)
  }, data = sim_data, method = method)

  rho_lower <- matrix(0, d, d)
  diag(rho_lower) <- 0.5
  rho_lower[lower.tri(rho_lower)] <- rho_lower_values
  rho_lower <- rho_lower + t(rho_lower)

  list(lower = rho_lower, upper = rho_upper)
}
