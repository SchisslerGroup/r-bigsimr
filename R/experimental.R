#' @export
convertCor_experimental <- function(rho,
                                    from = c("pearson", "spearman", "kendall"),
                                    to = c("pearson", "spearman", "kendall"),
                                    computeNearestPD = FALSE) {

  key <- paste(from, to)

  A <- switch (key,
               "pearson spearman" = function(r) (6/pi)*asin(r/2),
               "pearson kendall"  = function(r) (2/pi)*asin(r),
               "spearman pearson" = function(r) 2*sin(r*pi/6),
               "spearman kendall" = function(r) (2/pi)*asin(2*sin(r*pi/6)),
               "kendall pearson"  = function(r) sin(r*pi/2),
               "kendall spearman"  = function(r) (6/pi)*asin(sin(r*pi/2)/2),
               function(r) r
  )

  if (computeNearestPD) {
    return(Matrix::nearPD(A(rho), ensureSymmetry = FALSE, corr = TRUE)$mat)
  } else {
    return(A(rho))
  }
}


#' @export
cor2cor <- function(rho,
                    from = c("pearson", "spearman", "kendall"),
                    to = c("pearson", "spearman", "kendall")) {

  from <- match.arg(from)
  to   <- match.arg(to)

  CASE <- switch (paste(from, to, sep = "_"),
               "pearson_spearman" = 0L,
               "pearson_kendall"  = 1L,
               "spearman_pearson" = 2L,
               "spearman_kendall" = 3L,
               "kendall_pearson"  = 4L,
               "kendall_spearman" = 5L,
               -1L
  )

  CXX_cor2cor(rho, CASE)
}

#' @export
rvec_experimental <- function(n,
                              rho,
                              margins,
                              cores = 1,
                              type = c("pearson", "kendall", "spearman"),
                              computeNearestPD = TRUE,
                              useCholesky = FALSE,
                              useMod = FALSE) {

  # Internal function definitions
  .u2m <- function(u, margin) {
    margin$p <- quote(u)
    eval(margin)
  }

  .rmvuu <- function(n, rho, isChol = FALSE) {

    if (isChol) {

    }

    R    = reticulate::np_array(rho)

    R_d  = jax$device_put(R)
    m_d  = jax$numpy$zeros(reticulate::tuple(ncol(rho)))

    key  = jax$random$PRNGKey(sample(.Random.seed, 1))
    size = reticulate::tuple(as.integer(n))

    Z_d  = jax$random$multivariate_normal(key, m_d, R_d, size)$block_until_ready()
    U_d  = jax$scipy$stats$norm$cdf(Z_d)$block_until_ready()

    U    = jax$device_get(U_d)
    reticulate::py_to_r(U)

  }


}
