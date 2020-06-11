.u2m <- function(u, margin) {
  margin$p <- quote(u)
  eval(margin)
}
