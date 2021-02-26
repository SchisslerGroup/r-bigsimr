# Apply a set of stats to each column of a matrix
eachcol <- function(X, STATS, FUN) {
  sweep(X, 1, STATS, FUN)
}

# Apply a set of stats to each row of a matrix
eachrow <- function(X, STATS, FUN) {
  sweep(X, 2, STATS, FUN)
}
