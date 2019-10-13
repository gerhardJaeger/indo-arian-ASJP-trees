library(coda)

posterior1 <- read.table("output/indoiranian_run_1.log",
                         header=T,check.names=F)[,-1]
posterior2 <- read.table("output/indoiranian_run_2.log",
                         header=T,check.names=F)[,-1]
#posterior3 <- read.table("output/indoiranian_run_3.log",
#                         header=T,check.names=F)[,-1]
#posterior4 <-read.table("output/indoiranian_run_4.log",
#                        header=T,check.names=F)[,-1]

n <- min(nrow(posterior1),
         nrow(posterior2))

cn <- colnames(posterior1)

vcn <- cn[apply(rbind(posterior1[-(1:(n %/% 2)),],
                      posterior2[-(1:(n %/% 2)),]),
                2, var) > 0]


chains <- mcmc.list(list(mcmc(data=posterior1[1:n,vcn]),
                         mcmc(data=posterior2[1:n,vcn])))
gelman.diag(chains, multivariate = F)
