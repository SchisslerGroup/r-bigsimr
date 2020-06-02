#' Install bigsimr and its dependencies
#'
#' @inheritParams reticulate::conda_list
#'
#' @param method Installation method. By default, "auto" automatically finds a
#'   method that will work in the local environment. Change the default to force
#'   a specific installation method. Note that the "virtualenv" method is not
#'   available on Windows. Note also that since this command runs without
#'   privilege the "system" method is available only on Windows.
#'
#'
#' @param envname Name of Python environment to install within.
#'   Defaults to 'r-reticulate'
#'
#'
#' @param restart_session Restart R session after installing (note this will
#'   only occur within RStudio).
#'
#' @param conda_python_version the python version installed in the created conda
#'   environment. Python 3.6 is installed by default.
#'
#' @param ... other arguments passed to [reticulate::conda_install()] or
#'   [reticulate::virtualenv_install()].
#'
#' @export
install_bigsimr <- function(method = c("auto", "virtualenv", "conda"),
                            conda = "auto",
                            envname = NULL,
                            restart_session = TRUE,
                            conda_python_version = "3.7",
                            ...) {
  packages = c("jax", "jaxlib")

  method <- match.arg(method)

  reticulate::py_install(
    packages       = packages,
    envname        = envname,
    method         = method,
    conda          = conda,
    python_version = conda_python_version,
    pip            = TRUE,
    ...
  )

  cat("\nInstallation complete.\n\n")

  if (restart_session && rstudioapi::hasFun("restartSession"))
    rstudioapi::restartSession()

  invisible(NULL)
}
