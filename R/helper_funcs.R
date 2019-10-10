#' Converts Pearson correlation to Spearman. Used for normal random variables.
#'
#' @param R A a square symmetric Pearson correlation matrix.
#' @return A Spearman correlation matrix.
convertPearsonSpearman <- function(R) {
  tmp <- (6 / pi) * asin(R / 2)
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(R)
  return(tmp)
}


#' Converts Spearman correlation to Pearson. Used for normal random variables.
#'
#' @param R A a square symmetric Spearman correlation matrix.
#' @return A Pearson correlation matrix.
convertSpearmanPearson <- function(Rho) {
  tmp <- 2* sin( Rho * (pi / 6))
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(Rho)
  return(tmp)
}


#' Converts Pearson correlation to Kendall. Used for normal random variables.
#'
#' @param R A a square symmetric Pearson correlation matrix.
#' @return A Kendall correlation matrix.
convertPearsonKendall <- function(R) {
  tmp <- (2 / pi) * asin(R)
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(R)
  return(tmp)
}


#' Converts Kendall correlation to Pearson. Used for normal random variables.
#'
#' @param R A a square symmetric Kendall correlation matrix.
#' @return A Pearson correlation matrix.
convertKendallPearson <- function(Rho) {
  tmp <- sin( Rho * (pi / 2) )
  tmp <- as.matrix(Matrix::nearPD(tmp)$mat)
  rownames(tmp) <- colnames(tmp) <- colnames(Rho)
  return(tmp)
}


#' Adjust an input correlation matrix for our method to match
#'
#' @param R A a square symmetric correlation matrix.
#' @param params The parameters of the marginals.
#' @param cores The number of cores to utilize.
adjustInputR <- function(R, params, cores = 1){

  ## 1. find the pairs
  d <- NROW(R)
  index_mat <- combn(x = 1:d, 2)

  ## 2. compute the first order approximation
  if (cores == 1) {
    input_rho_vector <- apply(index_mat, 2, function(tmp_index){
      ## extract and compute simple quantities
      alpha1 <- params["alpha", tmp_index[1]]
      alpha2 <- params["alpha", tmp_index[2]]
      lambda_nb1 <- params["lambda", tmp_index[1]]
      lambda_nb2 <- params["lambda", tmp_index[2]]
      sd_nb1 <- sqrt(params["nb_var", tmp_index[1]])
      sd_nb2 <- sqrt(params["nb_var", tmp_index[2]])
      ## precompute scale factor
      scale_factor <- sqrt(alpha1) * sqrt(alpha2)
      ## extract target nb corr
      tmp_nb_rho <- R[tmp_index[1], tmp_index[2]]
      ## convert to rho for gamma (exact)
      target_gamma_rho <- tmp_nb_rho * ( (sd_nb1 * sd_nb2) / (scale_factor * lambda_nb1 * lambda_nb2) )
      ## now approximate rho for input corr for MVN
      return(min(1, target_gamma_rho))
    })
  } else {
    cl <- parallel::makeCluster(cores, type = "FORK")

    ## parallel version of the above code
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    input_rho_vector <- foreach::foreach(i = 1:ncol(index_mat), .combine = 'c') %dopar% {
      ## extract and compute simple quantities
      tmp_index <- index_mat[,i]
      ## extract and compute simple quantities
      alpha1 <- params["alpha", tmp_index[1]]
      alpha2 <- params["alpha", tmp_index[2]]
      lambda_nb1 <- params["lambda", tmp_index[1]]
      lambda_nb2 <- params["lambda", tmp_index[2]]
      sd_nb1 <- sqrt(params["nb_var", tmp_index[1]])
      sd_nb2 <- sqrt(params["nb_var", tmp_index[2]])
      ## precompute scale factor
      scale_factor <- sqrt(alpha1) * sqrt(alpha2)
      ## extract target nb corr
      tmp_nb_rho <- R[tmp_index[1], tmp_index[2]]
      ## convert to rho for gamma (exact)
      target_gamma_rho <- tmp_nb_rho * ((sd_nb1 * sd_nb2) / (scale_factor * lambda_nb1 * lambda_nb2))

      ## now approximate rho for input corr for MVN
      return(min(1, target_gamma_rho))
    }
    ## close cluster
    parallel::stopCluster(cl)
  }

  ## Put input_rho_vector into a matrix
  to_return <- diag(d)
  to_return[lower.tri(to_return)] <- input_rho_vector
  to_return <- to_return + t(to_return) - diag(diag(to_return))

  ## Ensure positive semi-definiteness
  eigen_obj <- eigen(to_return)
  if (any(eigen_obj$values < 0)) {
    to_return <- as.matrix(Matrix::nearPD(to_return)$mat)
  }

  ## Label and return
  rownames(to_return) <- colnames(to_return) <- rownames(R)
  return(to_return)
}


#' Transforms a [multivariate]normal vector to a different marginal via a
#' uniform intermediate transformation.
#'
#' @param x A normal random vector.
#' @param param A list ccontaining the marginal and its parameters.
normal2marginal <- function(x, param) {
  do.call(what = paste0("q", param[[1]]),
          args = c(list(p = pnorm(x)), param[-1]))
}
