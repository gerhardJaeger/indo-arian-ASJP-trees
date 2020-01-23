require(ggplot2)
require(ggrepel)
require(dplyr)
require(tikzDevice)
require(phangorn)
require(rstan)
require(reshape2)
require(plotluck)
require(Cairo)

orange <- "#E69F00"
pink <- "#CC79A7"
yellow <- "#F0E442"
red <- "#D55E00"
cyan <- "#56B4E9"
green <- "#009E73"
blue <- "#0072B2"

node.colors <- c(pink,orange,red,cyan,green,blue)

coords <- read.csv('language_locations.csv',header=F,sep=' ')
colnames(coords) <- c('lang','lon','lat')
coords <- coords[order(coords$lang),]

states <- read.csv('states.csv',header=F,sep=' ')
colnames(states) <- c('lang','state')
states <- states[order(states$lang),]

patmat <- read.csv('states_matter_pattern.csv',header=FALSE,sep='\t')
state.pat.mat<-paste(patmat$V1,patmat$V2,sep='')
names(state.pat.mat) <- patmat$V1

full.df <- cbind(coords,states$state)
colnames(full.df)[4] <- 'state'
full.df$lang.plain <- unlist(strsplit(as.character(full.df$lang),'<'))[seq(1,nrow(full.df)*2,2)]

clf <- unlist(strsplit(as.character(full.df$state),','))[1+c(1:nrow(full.df)-1)*2]
pl <- unlist(strsplit(as.character(full.df$state),','))[c(1:nrow(full.df))*2]

clf <- substr(clf,2,nchar(clf))
pl <- substr(pl,1,nchar(pl)-1)

pl[pl=='opt.pl'] <- '-obl.pl'
pl[pl=='obl.pl'] <- '+obl.pl'

full.df$pl <- as.factor(pl)
full.df$clf <- as.factor(clf)

full.df$state <- factor(full.df$state,levels=c('(-clf,morph.pl)','(-clf,opt.pl)','(-clf,obl.pl)','(+clf,opt.pl)','(+clf,obl.pl)'))

cairo_pdf('lang_map.pdf',height=5,width=10)
ggplot()+
borders(colour='lightgray',fill='lightgray')+coord_cartesian(xlim=c(-1.5,1.5)+range(coords$lon),ylim=c(-1.5,1.5)+range(coords$lat)) +

geom_point(data=full.df,aes(x=lon,y=lat,fill=state,color=state)) + scale_fill_manual(values=node.colors) + scale_color_manual(values=node.colors) +


#theme(legend.title = element_text(size = 14) ,legend.text = element_text(size = 14),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 14)) + 

#geom_point(data=full.df,aes(x=lon,y=lat,shape=state,fill=state,color=state)) + scale_fill_manual(values=node.colors) + scale_color_manual(values=node.colors) + scale_shape_manual(values = c(21:25)) +
#geom_text_repel(data=full.df,aes(x=lon,y=lat,label=lang.plain,color=state),cex=2)+scale_color_manual(values=node.colors)+
#geom_point(data=full.df,aes(x=lon,y=lat,fill=pl,color=pl,shape=clf)) + scale_fill_manual(values=node.colors) + scale_color_manual(values=node.colors) + scale_shape_manual(values = c(21:25)) +
#geom_text_repel(data=full.df,aes(x=lon,y=lat,label=lang.plain),cex=2)+
#geom_text_repel(data=full.df,aes(x=lon,y=lat,label=lang.plain,color=pl),cex=2)+scale_color_manual(values=node.colors)+
geom_text_repel(data=full.df,aes(x=lon,y=lat,label=lang.plain,color=state),cex=3)+scale_color_manual(values=node.colors)+
theme_bw()


dev.off()

#densitree

trees <- read.nexus('IIr.trees')

#cairo_pdf('tree_sample.pdf',height=15,width=10)

#densiTree(trees,alpha=.01)

#dev.off()

#rate graphs

load('all_rates.Rdata')

rates <- data.frame(extract(fit.full)$R)

state.values <- levels(states$V2)

state.combs <- c()
for (i in 1:length(state.values)) {
    for (j in 1:length(state.values)) {
        if (state.values[i]!=state.values[j]) {
            state.combs <- c(state.combs,paste(state.values[i],state.values[j],sep=' -> '))
        }
    }
}


state.combs <- state.combs[c(1,2,4,6,7,9,11,12)]

colnames(rates) <- state.combs

long.posterior <- melt(data=rates,variable.name='change.type',value.name='rate')

pdf('all_rates.pdf',colormodel='gray')
#plotluck(long.posterior,change.type~rate,color='black',fill='gray')
#ggplot(long.posterior, aes(x=rate, color=change.type)) + geom_density()
#ggplot(long.posterior, aes(x=rate, y=change.type)) + geom_density_ridges2()
#ggplot(long.posterior,aes(change.type,rate)) + geom_violin() + coord_flip()
ggplot(long.posterior,aes(change.type,rate)) + geom_violin(scale='width') + coord_flip() + labs(x='change type')
dev.off()

#diachronic GSS

dia.GSS <- data.frame(
    rates[,c('(-clf,-obl.pl) -> (+clf,-obl.pl)',
             '(-clf,+obl.pl) -> (+clf,+obl.pl)')]
)

colnames(dia.GSS) <- c('(-clf,-obl.pl) -> (+clf,-obl.pl)',
                       '(-clf,+obl.pl) -> (+clf,+obl.pl)')

#synchronic GSS
#exit rates

syn.GSS <- data.frame(
    rowSums(rates[,c(
        '(+clf,-obl.pl) -> (-clf,-obl.pl)',
        #'(+clf,-obl.pl) -> (-clf,+obl.pl)',
        '(+clf,-obl.pl) -> (+clf,+obl.pl)'
    )]),
    rowSums(rates[,c(
        #'(+clf,+obl.pl) -> (-clf,-obl.pl)',
        '(+clf,+obl.pl) -> (-clf,+obl.pl)',
        '(+clf,+obl.pl) -> (+clf,-obl.pl)'
    )])
)


colnames(syn.GSS) <- c('(+clf,-obl.pl)','(+clf,+obl.pl)')

#graphics

plot_multi_histogram <- function(df, feature, label_column, legend_label) {
    ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
        scale_fill_manual(values=c("#56B4E9","#CC79A7")) + 
        geom_histogram(alpha=0.7, position="identity", bins=150) +
        #geom_density(alpha=0.7) +
        geom_vline(aes(xintercept=median(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
        labs(x='rate', y = "") + 
        theme(axis.title.x = element_text(size = 22),axis.text.x = element_text(size = 22),axis.title.y = element_text(size = 22),legend.text=element_text(size=16)) +
        labs(fill=legend_label)	}


make_step <- function(df) {
    df_area_step <- bind_rows(old = df, new = df %>% mutate(counts = lag(counts)), .id = "source") %>% arrange(breaks, source)
    df_area_step[1,'counts'] <- 0
    return(df_area_step)
}


#raw distributions

syn.long <- melt(data=syn.GSS, variable.name='change.type',value.name='rate')
#syn.diff <- data.frame(difference.in.rates=syn.GSS[,'(+clf,+obl.pl)'] - syn.GSS[,'(+clf,-obl.pl)'])
#syn.diff <- syn.GSS[,'(+clf,+obl.pl)'] - syn.GSS[,'(+clf,-obl.pl)']
#syn.diff.hist <- data.frame(breaks=hist(syn.diff)$breaks,counts=c(hist(syn.diff)$counts,0))
#syn.diff.step <- make_step(syn.diff.hist)
diff.rates <- syn.GSS[,'(+clf,+obl.pl)'] - syn.GSS[,'(+clf,-obl.pl)']
print(paste(length(diff.rates[diff.rates<=0])/length(diff.rates),length(diff.rates[diff.rates>0])/length(diff.rates),length(diff.rates[diff.rates>0])/length(diff.rates[diff.rates<=0]),sep='|'))
diff.df <- with(density(diff.rates), data.frame(x,y))
colnames(diff.df) <- c('difference.in.rates','density')
N = length(diff.rates)
diff.rates <- sort(diff.rates)


cairo_pdf('syn_rates.pdf')
#plot_multi_histogram(syn.long[sample(c(1:nrow(syn.long)),10000),],"rate","change.type","exit rate")
plot_multi_histogram(syn.long,"rate","change.type","exit rate")
dev.off()

cairo_pdf('syn_diff.pdf')
ggplot(data = diff.df, mapping = aes(x = difference.in.rates, y = density)) + geom_line() +
    geom_area(mapping = aes(x = ifelse(difference.in.rates>0 ,difference.in.rates,0)), fill = "#56B4E9") + geom_area(mapping = aes(x = ifelse(difference.in.rates<=0 ,difference.in.rates,0)), fill = "#CC79A7") +
    labs(x='difference in rates') +
    geom_vline(aes(xintercept=0)) +
    theme(axis.title.x = element_text(size = 22),axis.text.x = element_text(size = 22),axis.title.y = element_text(size = 22)) +
    ylim(0,2) + xlim(-1,3)
dev.off()



dia.long <- melt(data=dia.GSS, variable.name='change.type',value.name='rate')
diff.rates <- dia.GSS[,'(-clf,-obl.pl) -> (+clf,-obl.pl)'] - dia.GSS[,'(-clf,+obl.pl) -> (+clf,+obl.pl)']
print(paste(length(diff.rates[diff.rates<=0])/length(diff.rates),length(diff.rates[diff.rates>0])/length(diff.rates),length(diff.rates[diff.rates>0])/length(diff.rates[diff.rates<=0]),sep='|'))
diff.df <- with(density(diff.rates), data.frame(x,y))
colnames(diff.df) <- c('difference.in.rates','density')
N = length(diff.rates)
diff.rates <- sort(diff.rates)


cairo_pdf('dia_rates.pdf')
#plot_multi_histogram(dia.long[sample(c(1:nrow(dia.long)),10000),],"rate","change.type","change type")
plot_multi_histogram(dia.long,"rate","change.type","change type")
dev.off()

cairo_pdf('dia_diff.pdf')
ggplot(data = diff.df, mapping = aes(x = difference.in.rates, y = density)) + geom_line() +
    geom_area(mapping = aes(x = ifelse(difference.in.rates>0 ,difference.in.rates,0)), fill = "#56B4E9") + geom_area(mapping = aes(x = ifelse(difference.in.rates<=0 ,difference.in.rates,0)), fill = "#CC79A7") +
    labs(x='difference in rates') +
    geom_vline(aes(xintercept=0)) +
    theme(axis.title.x = element_text(size = 22),axis.text.x = element_text(size = 22),axis.title.y = element_text(size = 22)) +
    ylim(0,2) + xlim(-1,3)
dev.off()