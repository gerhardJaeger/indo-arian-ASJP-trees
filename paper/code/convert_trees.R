require(phytools)

trees <- read.newick('indoiranian.posterior.tree')

states <- read.csv('states.csv',header=F,sep=' ',row.names=1)

key <- unlist(strsplit(rownames(states),split='<|>',perl=T))
key <- key[seq(2,length(key),2)]

states$value <- rownames(states)
rownames(states) <- key

new.trees <- list()

for (i in 1:length(trees)) {
	tree <- trees[[i]]
	for (j in 1:length(tree$node.label)) {
		if (tree$node.label[j] != '') {
			tree <- bind.tip(tree,tree$node.label[j],edge.length=1,where=length(tree$tip.label)+j)
		}
	}
    tree$tip.label <- states[tree$tip.label,]$value
	new.trees[[i]] <- tree
}

class(new.trees) <- 'multiPhylo'

write.nexus(file='IIr.trees',new.trees[seq(1,length(new.trees),10)])
