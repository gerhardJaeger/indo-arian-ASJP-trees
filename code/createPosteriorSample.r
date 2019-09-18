library(ape)

bi <- 0.5


trees1 <- read.nexus(paste0('output/indoiranian.run1.t'))
trees2 <- read.nexus(paste0('output/indoiranian.run2.t'))
trees3 <- read.nexus(paste0('output/indoiranian.run3.t'))
trees4 <- read.nexus(paste0('output/indoiranian.run4.t'))

trees <- c(trees1[-(1:floor(bi*length(trees1)))],
           trees2[-(1:floor(bi*length(trees2)))],
           trees3[-(1:floor(bi*length(trees3)))],
           trees4[-(1:floor(bi*length(trees4)))])
if (length(trees)>1000) {
    trees <- sample(trees,1000)
}
write.tree(trees, 'indoiranian.posterior.tree'))

