#' Computes the theoretical upper and lower bounds of possible correlations
#'   given a set of marginals
#'
#' @param margins The parameters of the marginals.
#' @param type The type of correlation matrix that is being passed.
#' @param cores The number of cores to utilize.
#' @param reps The number of reps used for computing the cor bounds
#' @return A list containing the theoretical upper and lower bounds
#' @export
computeCorBounds <- function(margins,
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
  sim_data <- sapply(1:d, function(i) {
    # Replace the quantile function with the RNG function (e.g. qnorm -> rnorm)
    margins[[i]][[1]] <- q2r(margins[[i]][[1]])
    margins[[i]]$n <- quote(reps)
    eval(rlang::call2("sort", margins[[i]]))
  })

  # Upper bounds
  rho_upper <- fastCor(sim_data, method = type)

  # Lower bounds
  rho_lower_values <-
    apply(index_mat, 2, function(index, data, ...) {
      fastCor(data[, index[1]],
              rev(data[, index[2]]),
              method = type)
    }, data = sim_data, method = type)

  rho_lower <- matrix(0, d, d)
  diag(rho_lower) <- 0.5
  rho_lower[lower.tri(rho_lower)] <- rho_lower_values
  rho_lower <- rho_lower + t(rho_lower)

  list(upper = rho_upper,
       lower = rho_lower)
}
