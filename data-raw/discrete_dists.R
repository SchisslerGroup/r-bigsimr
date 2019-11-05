## code to prepare `discrete_dists` dataset goes here
discrete_dists <- c(
  "binom",
  "geom",
  "hyper",
  "nbinom",
  "pois"
)

usethis::use_data(discrete_dists, internal = TRUE)
