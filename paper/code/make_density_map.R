require(rstan)
require(phytools)
require(phangorn)
require(expm)
require(tikzDevice)

IIr.trees <- read.nexus('IIr.trees')

#rate graphs

load('all_rates.Rdata')

mcc.tree <- maxCladeCred(IIr.trees)

q.obl.clf <- extract(fit.full)$R[,2]+extract(fit.full)$R[,3]
q.opt.clf <- extract(fit.full)$R[,5]+extract(fit.full)$R[,6]

all.rates <- data.frame(rate=c(q.obl.clf,q.opt.clf),change=c(rep('q.obl.clf',length(q.obl.clf)),rep('q.opt.clf',length(q.obl.clf))))

rates <- data.frame(extract(fit.full)$Q)
tree.ind <- rep(1:100,each=nrow(rates)/100)


orange <- "#E69F00"
pink <- "#CC79A7"
yellow <- "#F0E442"
red <- "#D55E00"
cyan <- "#56B4E9"
green <- "#009E73"
blue <- "#0072B2"

#node.colors <- c(pink,orange,red,cyan,green,blue)[2:5]
node.colors <- c(orange,red,cyan,green,blue)


F <- 4

fill.matrix <- function(R) {
    return(matrix(R,nrow=F,ncol=F))
}


states_n <- as.character(states[mcc.tree$tip.label,])
names(states_n) <- mcc.tree$tip.label

mcc.tree$edge.length <- mcc.tree$edge.length
mapped.trees <- list()
for (i in 1:1000) {
    print(i)
    t = sample(c(1:dim(rates)[1]),1)
    Q <- extract(fit.full)$Q[t,,]/1000
    rownames(Q) <- colnames(Q) <- levels(states$V2)
    mapped.trees[[i]] <- make.simmap(mcc.tree,states_n,Q=Q)
}

class(mapped.trees) <- c("multiSimmap","multiPhylo")

densiMap <- function(simmap.trees, pal=rainbow, alpha=10, ...){
# Palette is the color palette to draw from, alpha is the transparency of each plot. ... passed to plotSimmap.
states <- sort(unique(names(unlist(simmap.trees[[1]]$maps))))
nstates <- length(states)
makeTransparent <- function (someColor, alpha = 10) {
newColor <- col2rgb(someColor)
apply(newColor, 2, function(curcoldata) {
rgb(red = curcoldata[1], green = curcoldata[2], blue =
curcoldata[3],
alpha = alpha, maxColorValue = 255)
})
}
palC <- setNames(pal(nstates), states)
plotSimmap(simmap.trees[[1]], colors=makeTransparent(palC,alpha), ...)
dum <- lapply(2:length(simmap.trees), function(x)
plotSimmap(simmap.trees[[x]], colors=makeTransparent(palC,alpha),
add=TRUE,...))
}


densiMap <- function(simmap.trees, pal=node.colors, alpha=10, ...){
# Palette is the color palette to draw from, alpha is the transparency of each plot. ... passed to plotSimmap.
states <- sort(unique(names(unlist(simmap.trees[[1]]$maps))))
nstates <- length(states)
makeTransparent <- function (someColor, alpha = 10) {
newColor <- col2rgb(someColor)
apply(newColor, 2, function(curcoldata) {
rgb(red = curcoldata[1], green = curcoldata[2], blue =
curcoldata[3],
alpha = alpha, maxColorValue = 255)
})
}
palC <- setNames(pal, states)
plotSimmap(simmap.trees[[1]], colors=makeTransparent(palC,alpha), ...)
dum <- lapply(2:length(simmap.trees), function(x)
plotSimmap(simmap.trees[[x]], colors=makeTransparent(palC,alpha),
add=TRUE, ...))
}


#node.cols <- setNames(node.colors,sort(unique(states_n)))


#densiMap(mapped.trees, ftype="off", fsize=0.25)


boxlabel<-function(x,y,text,cex=1,bg="transparent",offset=0){
    w<-strwidth(text)*cex*1.1
    h<-strheight(text)*cex*1.4
    os<-offset*strwidth("W")*cex
    rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
    text(x,y,text,cex=cex,pos=4,offset=offset,font=1)
}


pdf('density_map.pdf',height=15,width=10)

par(fg="transparent")
densiMap(mapped.trees,
#ftype="off",
fsize=1.5, lwd=8)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
N<-Ntip(mcc.tree)
par(fg="black")
for(i in 1:Ntip(mcc.tree)) boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg=node.colors[as.numeric(states[mcc.tree$tip.label[i],])])


legend('topleft',legend=levels(states$V2),fill=node.colors,bty='n',cex=1.25)
dev.off()
