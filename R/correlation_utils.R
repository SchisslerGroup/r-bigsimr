#' Convert between different types of correlation matrices
#'
#' @param rho A a square symmetric correlation matrix
#' @param from The type of the input correlation matrix
#' @param to The type of the output correlation matrix
#' @export
convertCor <- function(rho,
                       from = c("pearson", "spearman", "kendall"),
                       to = c("pearson", "spearman", "kendall")) {

  key <- paste(from, to)

  A <- switch (key,
               "pearson spearman" = function(r) (6/pi)*asin(r/2),
               "pearson kendall"  = function(r) (2/pi)*asin(r),
               "spearman pearson" = function(r) 2*sin(r*pi/6),
               "spearman kendall" = function(r) (2/pi)*asin(2*sin(r*pi/6)),
               "kendall pearson"  = function(r) sin(r*pi/2),
               "kendal spearman"  = function(r) (6/pi)*asin(sin(r*pi/2)/2),
               function(r) r
  )

  Matrix::nearPD(A(rho), corr = TRUE)$mat
}


#' Computes the theoretical upper and lower bounds of possible correlations
#'   given a set of marginals
#'
#' @param margins The parameters of the marginals.
#' @param cores The number of cores to utilize.
#' @param type The type of correlation matrix that is being passed.
#' @return A list containing the theoretical upper and lower bounds
#' @export
computeCorBounds <- function(margins,
                             type = c("pearson", "kendall", "spearman"),
                             cores = 1,
                             reps = 1e5){

  type <- match.arg(type)
  d <- length(margins)
  index_mat <- combn(x = d, m = 2)


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
      eval( rlang::call2("sort", margins_new[[i]]) )
    })

    # Upper bounds
    rho_upper <- fastCor(sim_data, method = type)

    # Lower bounds
    rho_lower_values <- apply(index_mat, 2, function(index, data, ...){
      fastCor(data[, index[1]],
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
      eval( rlang::call2("sort", margins_new[[i]]) )
    }

    # Upper bounds (parallelized)
    rho_upper_values <- parallel::parApply(
      cl = cl,
      X = index_mat,
      MARGIN = 2,
      FUN = function(index, data, ...){
        fastCor(x = data[, index[1]],
                y = data[, index[2]],
                method = type)
      }, data = sim_data, method = type)

    rho_upper <- matrix(0, d, d)
    diag(rho_upper) <- 0.5
    rho_upper[lower.tri(rho_upper)] <- rho_upper_values
    rho_upper <- rho_upper + t(rho_upper)

    # Lower bounds
    rho_lower_values <- parallel::parApply(
      cl = cl,
      X = index_mat,
      MARGIN = 2,
      FUN = function(index, data, ...){
        fastCor(x = data[, index[1]],
                y = rev(data[, index[2]]),
                method=type)
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


#' Use the fastest methods available for computing each type of correlation
#'
#' @param x A matrix or vector
#' @param y A vector
#' @param method The type of correlation to compute
#' @export
fastCor <- function(x, y = NULL, method = c("pearson", "kendall", "spearman")) {

  method <- match.arg(method)

  stopifnot(method %in% c("pearson", "kendall", "spearman"))

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (method == "pearson") {
    if (is.null(y)) {
      coop::pcor(x)
    } else {
      coop::pcor(x, y)
    }

  } else if (method == "spearman") {
    if (is.null(y)) {
      coop::pcor(apply(x, 2, fastrank::fastrank_average))
    } else {
      coop::pcor(fastrank::fastrank_average(x),
                 fastrank::fastrank_average(y))
    }

  } else {
    if (is.null(y)) {
      pcaPP::cor.fk(x)
    } else {
      pcaPP::cor.fk(x, y)
    }

  }

}
