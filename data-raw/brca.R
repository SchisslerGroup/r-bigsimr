## code to prepare `brca` dataset goes here
library(tidyverse)
brca <- read.csv("data-raw/brca200.csv")
id <- brca[,1]
brca <- brca[,-1]
rownames(brca) <- id

usethis::use_data(brca, overwrite = TRUE)
