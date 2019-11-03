library(coda)

modelName <- "gamma11"

## posterior1 <- read.table(paste0("output/", modelName, "_run_1.p"),
##                          header=T,check.names=F)
## posterior2 <- read.table(paste0("output/", modelName, "_run_2.p"),
##                          header=T,check.names=F)
## posterior3 <- read.table(paste0("output/", modelName, "_run_3.p"),
##                          header=T,check.names=F)
## posterior4 <- read.table(paste0("output/", modelName, "_run_4.p"),
##                          header=T,check.names=F)

posterior <- read.table(paste0("output/", modelName, ".p"),
                        header=T, check.names=F)

posterior1 <- posterior[posterior$Replicate_ID==0, -2]
posterior2 <- posterior[posterior$Replicate_ID==1, -2]
posterior3 <- posterior[posterior$Replicate_ID==2, -2]
posterior4 <- posterior[posterior$Replicate_ID==3, -2]


n <- min(nrow(posterior1),
         nrow(posterior2),
         nrow(posterior3),
         nrow(posterior4))


chains <- mcmc.list(list(mcmc(data=posterior1[1:n,-1]),
                         mcmc(data=posterior2[1:n,-1]),
                         mcmc(data=posterior3[1:n,-1]),
                         mcmc(data=posterior4[1:n,-1])))


gelman.diag(chains, autoburnin=F, transform=T, multivariate=F)
