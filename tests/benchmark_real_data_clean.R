library(tidyverse)

# Load in the data
dat0 <- readRDS(
  file = "~/Downloads/complete_processed_tcga2stat_RNASeq2_with_clinical.rds")

# Select just the BRCA genes
brca0 <- dat0 %>%
  filter(disease == "BRCA") %>%
  select(-c("patient":"tumorsize"))

brca1 <- round(brca0, 0)

rm(dat0)

# Filter down the brca data set
top_percentile <- 0.95
brca_median <- apply(brca1, 2, median)
cut_point <- quantile(brca_median, top_percentile)
keep_genes <- names(brca_median)[brca_median >= cut_point]

(d <- length(keep_genes))

brca2 <- brca1 %>%
  select(all_of(keep_genes))


# Estimate the Spearman correlation
brca_rho <- bigsimr::fastCor(brca2, method = "spearman")

# Estimate the marginal parameters assuming negative binomial distribution
mom_nbinom <- function(x) {
  m <- mean(x)
  s <- sd(x)
  list("nbinom", size = m^2 / (s^2 - m), prob = m / s^2)
}

brca_margins <- apply(brca2, 2, mom_nbinom)

saveRDS(brca_rho, "tests/brca_rho.rds")
saveRDS(brca_margins, "tests/brca_margins.rds")
saveRDS(brca2, "tests/brca.rds")

brca_rho <- readRDS("tests/brca_rho.rds")
brca_margins <- readRDS("tests/brca_margins.rds")
brca2 <- readRDS("tests/brca.rds")

# Go through rvec algorithm --------------------------------------------------
devtools::load_all()
library(profvis)

n <- nrow(brca2)


p <- profvis({
  x <- rvec(n, brca_rho, brca_margins,
            cores = 1,
            type = "spearman")
})
