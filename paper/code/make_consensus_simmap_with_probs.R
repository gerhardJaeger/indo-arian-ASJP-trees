require(rstan)
require(phytools)
require(phangorn)
require(expm)
#require(tikzDevice)

IIr.trees <- read.nexus('IIr.trees')

#rate graphs

load('all_rates.Rdata')

mcc.tree <- maxCladeCred(IIr.trees)

q.obl.clf <- extract(fit.full)$R[,2]+extract(fit.full)$R[,3]
q.opt.clf <- extract(fit.full)$R[,5]+extract(fit.full)$R[,6]

all.rates <- data.frame(rate=c(q.obl.clf,q.opt.clf),change=c(rep('q.obl.clf',length(q.obl.clf)),rep('q.opt.clf',length(q.obl.clf))))

rates <- data.frame(extract(fit.full)$R)
tree.ind <- rep(1:100,each=nrow(rates)/100)


orange <- "#E69F00"
pink <- "#CC79A7"
yellow <- "#F0E442"
red <- "#D55E00"
cyan <- "#56B4E9"
green <- "#009E73"
blue <- "#0072B2"

#node.colors <- c(pink,orange,red,cyan,green,blue)[2:5]
node.colors <- c(orange,red,cyan,green)


F <- 4

fill.matrix <- function(R) {
    Q <- matrix(nrow=F,ncol=F)
    k <- 1
    for (i in 1:F) {
        for (j in 1:F) {
            if (i != j) {
                Q[i,j] <- R[[k]]
                k <- k + 1
            }
        }
    }
    for (i in 1:F) {
        z <- 0
        for (j in 1:F) {
            if (i != j) {
                z <- z - Q[i,j]
            }
        }
        Q[i,i] <- z
    }
    return(Q)
}


states_n <- as.character(states[mcc.tree$tip.label,])
names(states_n) <- mcc.tree$tip.label


mapped.trees <- list()
for (i in 1:1000) {
	print(i)
	t = sample(c(1:dim(rates)[1]),1)
    Q <- extract(fit.full)$Q[t,,]/1000
	rownames(Q) <- colnames(Q) <- levels(states$V2)
	mapped.trees[[i]] <- make.simmap(mcc.tree,states_n,Q=Q)
}

class(mapped.trees) <- c("multiSimmap","multiPhylo")



transitions <- list()
for (j in 1:length(mapped.trees[[1]]$maps)) { transitions[[j]] <- list() }

for (i in 1:length(mapped.trees)) {
	for (j in 1:length(mapped.trees[[i]]$maps)) {
		transitions[[j]][[i]] <- paste(names(mapped.trees[[i]]$maps[[j]]),collapse=' ')
	}
}

MAP <- list()
proportion <- list()
for (j in 1:length(mapped.trees[[1]]$maps)) {
	MAP[[j]] <- names(sort(table(unlist(transitions[[j]])),decreasing=TRUE))[1]
	proportion[[j]] <- sort(table(unlist(transitions[[j]])),decreasing=TRUE)[1]
}


MAP.map <- list()
for (i in 1:length(MAP)) {
	states.i <- unlist(strsplit(MAP[[i]],' '))
	n_states <- length(states.i)
	chunks <- mcc.tree$edge.length[[i]]/n_states
	intervals <- rep(chunks,n_states)
	names(intervals) <- states.i
	MAP.map[[i]] <- intervals
}


boxlabel<-function(x,y,text,cex=1,bg="transparent",offset=0){
    w<-strwidth(text)*cex*1.1
    h<-strheight(text)*cex*1.4
    os<-offset*strwidth("W")*cex
    rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
    text(x,y,text,cex=cex,pos=4,offset=offset,font=1)
}




final.map <- mcc.tree

class(final.map) <- c('simmap','phylo')

final.map$maps <- MAP.map

node.cols <- setNames(node.colors,sort(unique(states_n)))

pdf('consensus_simmap.pdf',height=15,width=10)
par(fg="transparent")


plot(final.map,col=node.cols,lwd=8)

pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
N<-Ntip(mcc.tree)
par(fg="black")
for(i in 1:Ntip(mcc.tree)) boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg=node.colors[as.numeric(states[mcc.tree$tip.label[i],])])
edgelabels(unlist(proportion)/1000,bg=NULL,frame='none',cex=.75,adj = c(0.5, -0.25))


legend('topleft',legend=levels(states$V2),fill=node.colors,bty='n',cex=1.25)
dev.off()
