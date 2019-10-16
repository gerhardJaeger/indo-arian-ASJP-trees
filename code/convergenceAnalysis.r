library(ape)
library(phytools)
library(coda)

p1 <- read.table('output/indoiranian_run_1.trees',
                 header=T, as.is=T)
#l1 <- nrow(p1)
#p1 <- p1[(l1-3999):l1,]

p2 <- read.table('output/indoiranian_run_2.trees',
                 header=T, as.is=T)
#l2 <- nrow(p2)
#p2 <- p2[(l2-3999):l2,]


tree = read.tree(text=p1$tree[1])

taxa <- sort(tree$tip.label)

    
trees1 <- lapply(p1$tree,
                 function(x)
                     rotateConstr(read.tree(text=x),
                                  taxa))
trees2 <- lapply(p2$tree,
                 function(x)
                     rotateConstr(read.tree(text=x),
                                  taxa))


getNodeHeights <- function(tree) {
    eHeights <- nodeHeights(tree)
    n <- max(tree$edge)
    depths <- sapply(1:n,
                     function(x)
                         eHeights[tree$edge == x][1])
    output <- rep(0, 127)
    output[1:length(depths)] <- max(depths)-depths
    return(output)
}

c1 <- t(sapply(trees1, getNodeHeights))
c2 <- t(sapply(trees2, getNodeHeights))
colnames(c1) <- paste('node', 1:127, sep='')
colnames(c2) <- paste('node', 1:127, sep='')


n <- min(nrow(c1), nrow(c2))

chains <- mcmc.list(list(data=mcmc(c1[1:n,]),
                         data=mcmc(c2[1:n,])))

gelman.diag(chains, multivariate = T, autoburnin = F)

gelman.plot(chains, autoburnin = F)


########

getClades <- function(tree) {
    st <- subtrees(tree)
    clades <- sapply(st, function(x) x$tip.label)
    splits <- unique(sapply(clades, function(x) paste0(1*(taxa %in% x), collapse='')))
    return(splits)
}

cf1 <- c()
for (tree in trees1[1:n]) {
    splits <- getClades(tree)
    for (s in splits) {
        if (s %in% names(cf1)) {
            cf1[s] <- cf1[s] + 1
        } else {
            cf1[s] <- 1
        }
    }
}

cf2 <- c()
for (tree in trees2[1:n]) {
    splits <- getClades(tree)
    for (s in splits) {
        if (s %in% names(cf2)) {
            cf2[s] <- cf2[s] + 1
        } else {
            cf2[s] <- 1
        }
    }
}

minpartfreq <- 0.1

partitions <- unique(c(names(cf1), names(cf2)))

for (p in partitions) {
    if (!(p %in% names(cf1))) {
        cf1[p] <- 0
    }
    if (!(p %in% names(cf2))) {
        cf2[p] <- 0
    }
}


partition.frequencies <- data.frame(chain1 = cf1[partitions],
                                    chains2 = cf2[partitions])

partition.frequencies$sum <- apply(partition.frequencies, 1, sum) / (2*n)

partition.frequencies <- partition.frequencies[partition.frequencies$sum > minpartfreq,]
rownames(partition.frequencies) <- 1:nrow(partition.frequencies)


print(mean(apply(partition.frequencies[,1:2] / n, 1, sd)))

for (t in c(trees1[1:n], trees2[1:n])) {
    write.tree(t, file = 'indoiranian.posterior.tree', append=T)
}
