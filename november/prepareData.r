library(ape)

d <- read.table("states.csv", as.is=T)
colnames(d) <- c('language', 'state')

## 1: (+clf,obl.pl)
## 2: (+clf,opt.pl)
## 3: (-clf,obl.pl)
## 4: (-clf,opt.pl)
## 5: (-clf,morph.pl)

## direct transitions:
## 1 -> 2
## 1 -> 3
## 2 -> 1
## 2 -> 4
## 3 -> 1
## 3 -> 4
## 3 -> 5
## 4 -> 2
## 4 -> 3
## 4 -> 5
## 5 -> 3
## 5 -> 4

## indirekt transitions:
## 1 -> 4
## 1 -> 5
## 2 -> 3
## 2 -> 5
## 3 -> 2
## 4 -> 1
## 5 -> 1
## 5 -> 2


key <- 1:5
names(key) <- c("(+clf,obl.pl)",
                "(+clf,opt.pl)",
                "(-clf,obl.pl)",
                "(-clf,opt.pl)",
                "(-clf,morph.pl)")

d$code <- key[d$state]

d$taxon <- sapply(strsplit(d$language, "[<>]"), function(x) x[2])



nexusHeader = paste0("#Nexus\n
Begin Data;\n
Dimensions ntax=", nrow(d),
" nchar=1;\n
Format Datatype=standard gap=? missing=- symbols=\"12345\";\n
Matrix")

nexusFooter = ";\nEnd;\n"

write(nexusHeader, 'states.nex')
write.table(d[,c('taxon', 'code')],
            'states.nex',
            quote=F,
            row.names=F,
            col.names=F,
            sep='\t',
            append=T)
write(nexusFooter, 'states.nex', append=T)


trees <- read.nexus('IIr.trees')



for (i in 1:length(trees)) {
    tr <- trees[[i]]
    tr$tip.label <- sapply(strsplit(tr$tip.label, "[<>]"), function(x) x[2])
    write.tree(tr, 'IIr.tree', append=T)
}

