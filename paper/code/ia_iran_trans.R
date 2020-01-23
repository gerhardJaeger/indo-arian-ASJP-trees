require(rstan)
require(phytools)
require(phangorn)
require(expm)
require(reshape2)
require(tikzDevice)

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
node.colors <- c(orange,red,cyan,green,blue)


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


ir.opt <- c()
ir.obl <- c()
ia.opt <- c()
ia.obl <- c()



get.anc.edges <- function(x,tip) {
    anc.nodes <- Ancestors(x,tip)
    anc.edges <- rbind(c(anc.nodes[1],tip),cbind(anc.nodes[2:length(anc.nodes)],anc.nodes[1:(length(anc.nodes)-1)]))
    #edge.inds <- intersect(x$edge,anc.edges)
    edge.inds <- match(paste(anc.edges[,1],anc.edges[,2]),paste(x$edge[,1],x$edge[,2]))
    return(edge.inds)
}


get.lang.history <- function(x,inds) {
    states.h <- c()
    for (i in 1:length(inds)) {
        states.h <- c(states.h,rev(names(x$maps[[inds[i]]])))
    }
    #print(paste(states.h,collapse=';'))
    prev.state <- ''
    for (i in 1:length(states.h)) {
        if (substr(states.h[i],1,2) == '(-') {
            prev.state <- states.h[i]
            break
        }
    }
    return(prev.state)
}


clf.langs <- rownames(states)[substr(as.character(states$V2),1,2)=='(+']
langs <- c()
prev.states <- c()

for (i in 1:1000) {
	print(i)
	tree <- sample(IIr.trees,1)[[1]]
    #ir.tree <- extract.clade(tree,100)
    #ia.tree <- extract.clade(tree,67)
	t = sample(c(1:dim(rates)[1]),1)
	#rate.slice <- as.vector(rates[t,])
	Q <- extract(fit.full)$Q[t,,]/1000
	#Q <- fill.matrix(rate.slice)/1000
	rownames(Q) <- colnames(Q) <- levels(states$V2)
    #ir.mapped.tree <- make.simmap(ir.tree,states_n[ir.tree$tip.label],Q=Q)
    #ia.mapped.tree <- make.simmap(ia.tree,states_n[ia.tree$tip.label],Q=Q)
	mapped.tree <- make.simmap(tree,states_n,Q=Q)
    ir.mapped.tree <- extract.clade.simmap(mapped.tree,100)
	ir.obl <- c(ir.obl,sum(summary(ir.mapped.tree)$Tr['(-clf,+obl.pl)',3:4]))
	ir.opt <- c(ir.opt,sum(summary(ir.mapped.tree)$Tr['(-clf,-obl.pl)',3:4]))
	ia.mapped.tree <- extract.clade(mapped.tree,67)
	ia.obl <- c(ia.obl,sum(summary(ia.mapped.tree)$Tr['(-clf,+obl.pl)',3:4]))
	ia.opt <- c(ia.opt,sum(summary(ia.mapped.tree)$Tr['(-clf,-obl.pl)',3:4]))
    for (j in 1:length(clf.langs)) {
        lang <- clf.langs[j]
        lang.edges <- get.anc.edges(mapped.tree,which(mapped.tree$tip.label==lang))
        prev.state <- get.lang.history(mapped.tree,lang.edges)
        langs <- c(langs,lang)
        prev.states <- c(prev.states,prev.state)
    }
    
}

lang.order <- clf.langs[order(match(clf.langs,mapped.tree$tip.label))]

prev.states <- data.frame('lang'=langs,'prev.state'=prev.states)
prev.state.probs <- as.data.frame.matrix(prop.table(xtabs( ~ lang + prev.state, prev.states),1))
prev.state.probs <- prev.state.probs[lang.order,c('(-clf,-obl.pl)','(-clf,+obl.pl)')]

patmat <- read.csv('states_matter_pattern.csv',header=FALSE,sep='\t')
state.pat.mat<-paste(patmat$V1,patmat$V2,sep='')
names(state.pat.mat) <- patmat$V1

new.row.names <- unname(state.pat.mat[rownames(prev.state.probs)])
new.row.names <- gsub('<.+>','',new.row.names,perl=TRUE)
new.row.names <- gsub('_',' ',new.row.names,perl=TRUE)
new.row.names <- gsub('\\*','$*$',new.row.names,perl=FALSE)
new.row.names <- gsub('†','$\\\\dagger$',new.row.names,perl=FALSE)

rownames(prev.state.probs) <- new.row.names


write.table(prev.state.probs,file='prev_state.csv',quote=FALSE,sep=' & ',eol='\\\\\n',row.names=TRUE,col.names=TRUE)


#ia.data <- as.data.frame(xtabs( ~ ia.opt + ia.obl))
#ggplot(ia.data,aes(ia.opt,ia.obl,fill=Freq))+geom_tile()


#ir.data <- as.data.frame(xtabs( ~ ir.opt + ir.obl))
#ggplot(ir.data,aes(ir.opt,ir.obl,fill=Freq))+geom_tile()

#p<- ggplot() + geom_tile(data=ir.data,aes(x=ir.opt,y=ir.obl,fill=Freq)) + scale_fill_gradient(low='white',high='blue') +
#    new_scale_fill() +
#    geom_tile(data=ia.data,aes(x=ia.opt,y=ia.obl,fill=Freq)) +
#    scale_fill_gradient(low='white',high='red') +
#    new_scale_fill()


#ggplot(ia.data,aes(x=ia.opt,y=ia.obl)) +
#theme_bw()+
#geom_tile(alpha=ir.data$Freq[1:247],fill='blue')+
#geom_tile(alpha=ia.data$Freq,fill='red')+
#theme(panel.grid.major=element_blank())



#ggplot(iran.long,aes(x=number.of.transitions)) + 
#    geom_histogram(data=subset(iran.long,change.type == "-clf -> +clf | -obl.pl"),fill = "red", alpha = 0.2) +
#    geom_histogram(data=subset(iran.long,change.type == "-clf -> +clf | +obl.pl"),fill = "blue", alpha = 0.2)




plot_multi_histogram <- function(df, feature, label_column) {
    ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
	scale_fill_manual(values=c("#56B4E9","#CC79A7")) + 
    geom_histogram(alpha=0.7, position="identity", bins=15) +
    #geom_density(alpha=0.7) +
    #geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x='number of transitions', y = "count") + 
	theme(axis.title.x = element_text(size = 22),axis.text.x = element_text(size = 22),axis.title.y = element_text(size = 22),legend.text=element_text(size=16)) +
	    labs(fill='change type')	}



ir.transitions <- data.frame(ir.opt,ir.obl)
colnames(ir.transitions) <- c("-clf -> +clf | -obl.pl","-clf -> +clf | +obl.pl")
ir.long <- melt(data=ir.transitions,variable.name='change.type',value.name='number.of.transitions')


cairo_pdf('iran_transitions.pdf')
plot_multi_histogram(ir.long, 'number.of.transitions', 'change.type')
dev.off()




ia.transitions <- data.frame(ia.opt,ia.obl)
colnames(ia.transitions) <- c("-clf -> +clf | -obl.pl","-clf -> +clf | +obl.pl")
ia.long <- melt(data=ia.transitions,variable.name='change.type',value.name='number.of.transitions')

cairo_pdf('ia_transitions.pdf')
plot_multi_histogram(ia.long, 'number.of.transitions', 'change.type')
dev.off()
