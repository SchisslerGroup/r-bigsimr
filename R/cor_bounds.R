#' Computes the theoretical upper and lower bounds of possible correlations
#'   given a set of marginals
#'
#' @param margins A list containing the marginal distributions which are
#'     distributions from the `distr6` package
#' @param method The type of correlation matrix that is being passed.
#' @param reps The number of sims used to estimate the bounds.
#' @return A list containing the theoretical upper and lower bounds
#' @export
cor_bounds <- function(margins,
                       method = c("pearson", "kendall", "spearman"),
                       reps = 1e5){

  method <- match.arg(method)

  # Because the margins are a list of `distr6` distributions (essentially envs)
  # and because of the S3 like nature of `rand`, a "new" rand function must
  # be made within this calling environment so that the number of reps gets
  # passed correctly (by baking it into the function).
  eval(substitute(local_rand <- function(X) distr6::rand(X, n = r),
                  list(r = reps)))

  # Generate random samples for each margin and sort the vectors
  # The sorted vectors are used for computing the bounds
  simdata <- sapply(margins, local_rand)
  simsort <- apply(simdata, 2, sort)
  simrev  <- apply(simsort, 2, rev)

  upper <- cor_fast(simsort, method = method)
  lower <- cor(simsort, simrev, method = method)
  diag(lower) <- 1.0

  list(lower = lower, upper = upper)
}
