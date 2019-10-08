source("R/crude_helper_funcs.R")
source("R/helper_funcs.R")
source("R/simInvNB.R")

simInvNB_old <- function(n, R, params, cores = 1, type = c("pearson", "kendall", "spearman")){
    
    ## Handle different types of dependencies
    if (type == "spearman") {
        R <- convertSpearmanPearson(R)
    }
    if (type == "kendall") {
        R <- convertKendallPearson(R)
    }
    
    ## determine the dimension d
    d <- NROW(R)

    ## 1. generate MVN sample
    mu <- rep(0, d)
    mvn_sim <- mvnfast::rmvn(n = n, mu = mu, sigma = R, ncores = cores, isChol = FALSE)
    ## cor(mvn_sim)

    ## single threaded
    if ( cores == 1 ) {
        ## 2. Apply normal cdf to obtain copula
        unif_sim <- apply(mvn_sim, 2, pnorm)
        remove(mvn_sim)
        ## cor(unif_sim)
        
        ## 3. now DIRECTLY TO negative binomial inverse transform to gamma
        nb_sim <- sapply(1:d, function(i){
            my_nb_prob <- (1 + params["lambda",i])^(-1)
            qnbinom(p = unif_sim[,i], size = params["alpha",i] , prob = my_nb_prob)
        })
        remove(unif_sim)
        ## cor(nb_sim); R

    } else {
        ## mulit-core --- parallelize
        ## 1. Start cluster
        cl <- parallel::makeCluster(cores, type = "FORK")

        ## 2. Apply normal cdf to obtain copula
        unif_sim <- parallel::parApply(cl = cl, X = mvn_sim, 2, pnorm)
        remove(mvn_sim)
        
        ## 3. now DIRECTLY TO  transform to negative binomial
        doParallel::registerDoParallel(cl)
        nb_sim <- foreach::foreach(i = 1:d, .combine = 'cbind') %dopar% {
            my_nb_prob <- (1 + params["lambda",i])^(-1)
            qnbinom(p = unif_sim[,i], size = params["alpha",i] , prob = my_nb_prob)
        }
        remove(unif_sim)

        ## 4. close cluster
        parallel::stopCluster(cl)
    }
    
    ## 5. Return the simulated data set
    colnames(nb_sim) <- rownames(R)
    return(nb_sim)
}

library(ggplot2)
library(microbenchmark)
library(foreach)

#= D=10, N=1e5, Cores=24
d <- 10
set.seed(12)
rho <- rcorr(d = d)
params <- rnegbin_params(d = d)
n <- 1e5
cores <- 10

mbm <- microbenchmark(
    "0simInvNB" = simInvNB_old(n, rho, params, 1, "pearson"),
    "1simInvNB_new" = simInvNB(n, rho, params, 1, "pearson"),
    "2simInvNB_24core" = simInvNB_old(n, rho, params, cores, "pearson"),
    "3simInvNB_new_24core" = simInvNB(n, rho, params, cores, "pearson")
)

(p <- autoplot(mbm) + labs(title = "dims=10, n=1e5, cores=1, 24"))
ggsave("benchmark_simInvNB_d10n1e5c24.png", p, device = "png")


#= D=100, N=1e6, Cores=24
d <- 100
set.seed(12)
rho <- rcorr(d = d)
params <- rnegbin_params(d = d)
n <- 1e6
cores <- 24

mbm <- microbenchmark(
    "0simInvNB" = simInvNB_old(n, rho, params, 1, "pearson"),
    "1simInvNB_new" = simInvNB(n, rho, params, 1, "pearson"),
    "2simInvNB_24core" = simInvNB_old(n, rho, params, cores, "pearson"),
    "3simInvNB_new_24core" = simInvNB(n, rho, params, cores, "pearson")
)
(p <- autoplot(mbm) + labs(title = "dims=100, n=1e6, cores=1, 24"))
ggsave("benchmark_simInvNB_d100n1e6c24.png", p, device = "png")