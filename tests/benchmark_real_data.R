library(profvis)
library(tidyverse)
# library(bigsimr)
devtools::load_all()
CORES <- parallel::detectCores()
set.seed(06082020)

## full workflow simulating RNA-seq data
allDat <- readRDS(
  file = "~/Downloads/complete_processed_tcga2stat_RNASeq2_with_clinical.rds")
lastClinical <- which( names(allDat) == 'tumorsize' )
brca0 <- allDat[ allDat$disease == "BRCA", (lastClinical+1):ncol(allDat) ]

## remove naively the RSEM adjustment
brca1 <- round(brca, 0)
ncol(brca1) ## num of genes 20501

## compute the average expression
brcaMedian <- apply(brca1, 2, median)
## retain top (1-probs)*100% highest expressing genes for illustration
myProb <- 0.95
cutPoint <- quantile( x = brcaMedian, probs = myProb )
genesToKeep <- names( brcaMedian ) [ which(brcaMedian >= cutPoint) ]
brca2 <- brca1[ , genesToKeep ]
## ncol(brca) / 20501
(d <- length(genesToKeep))


logVarOverMean <- function( x ) {
  log ( var(x) / mean (x) )
}

summaryBRCA <- brca2 %>%
  summarize_all(  list(~ logVarOverMean( . ) ) )

summaryBRCA <- as.data.frame( t(as.data.frame( summaryBRCA )) )

names(summaryBRCA) <- "logVarOverMean"

ggplot(data = summaryBRCA, mapping = aes(x = logVarOverMean)) +
  geom_histogram(color = "white", bins = 20)

corType <- 'spearman'
system.time(rho <- bigsimr::fastCor(brca, method = corType))

estimateNegBinMoM <- function(tmpGene, minP = 1e-5) {
  tmpMean <- mean(tmpGene)
  tmpVar <- var(tmpGene)
  ## relate to nbinom parameters
  ## See ?rbinom for details.
  p <- tmpMean / tmpVar
  if (p < minP) {p <- minP } ## maybe add noise here
  n <- ( tmpMean * p) / ( 1  - p )
  ## format for bigsimr margins
  return( list("nbinom", size = n, prob = p) )
}

brcaMargins <- apply( unname(as.matrix(brca)), 2, estimateNegBinMoM )

N <- nrow(brca)

system.time({
  simBRCA <- rvec(N,
                  rho = rho,
                  params = brcaMargins,
                  cores = 1,
                  type = corType,
                  adjustForDiscrete = FALSE)
})

system.time({
  simBRCA <- rvec(N,
                  rho = rho,
                  params = brcaMargins,
                  cores = CORES,
                  type = corType,
                  adjustForDiscrete = FALSE)
})

profvis({
  simBRCA <- rvec(N,
                  rho = rho,
                  params = brcaMargins,
                  cores = 1,
                  type = corType,
                  adjustForDiscrete = FALSE)
})

system.time({
  convertCor(rho, "spearman", "pearson")
})

