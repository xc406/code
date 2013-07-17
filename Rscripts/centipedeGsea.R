\section[examples]{Examples\footnote{examples' note}}

comparisons
YY1/GATA2 
1. peak-centric comparison
centipede to luciferase
htseq to luciferase
2. gene-centric comparison
sumOfReads to luciferase
poisson to luciferase


penalized for false negatives

sometimes, gs does not have enough coverage, true negatives might not be true

rm(list=ls())

source('~/Documents/gitCode/Rscripts/gsea.R')
source('~/Documents/gitCode/Rscripts/aupr.R')

file.yy1 <- "~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Yy1sc281V0416101PkRep1genes_format"
pred.gs.yy1 <- read.delim(file.yy1,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
pred.yy1 <- sort(pred.yy1,decreasing=T)
gs.yy1 <- gs.yy1>0
AUPR.res.yy1 <- calcAupr(pred.yy1,gs.yy1)

file.gata2 <- "~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Gata2sc267Pcr1xPkRep1genes_format"
pred.gs.gata2 <- read.delim(file.gata2,header=F)
pred.gata2 <- -log10(ppois(pred.gs.gata2[,3],pred.gs.gata2[,4],lower.tail=F))
gs.gata2 <- pred.gs.gata2[,2]
names(gs.gata2) <- pred.gs.gata2[,1]
names(pred.gata2) <- names(gs.gata2)
pred.gata2 <- sort(pred.gata2,decreasing=T)
gs.gata2 <- gs.gata2>0
AUPR.res.gata2 <- calcAupr(pred.gata2,gs.gata2)

file.max <- "~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562MaxV0416102PkRep1genes_format"
pred.gs.max <- read.delim(file.max,header=F)
pred.max <- -log10(ppois(pred.gs.max[,3],pred.gs.max[,4],lower.tail=F))
gs.max <- pred.gs.max[,2]
names(gs.max) <- pred.gs.max[,1]
names(pred.max) <- names(gs.max)
pred.max <- sort(pred.max,decreasing=T)
gs.max <- gs.max>0
AUPR.res.max <- calcAupr(pred.max,gs.max)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~comparing function

file.yy1 <- "~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_htseq_gname_format_cut"
pred.gs.yy1 <- read.delim(file.yy1,header=F)
pred.yy1 <- pred.gs.yy1[,2]
gs.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
#pred.yy1 <- sort(pred.yy1,decreasing=T)
gs.yy1 <- gs.yy1>0
AUPR.res.yy1 <- calcAupr(pred.yy1,gs.yy1)

file.gata2 <- "~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_htseq_gname_format_cut"
pred.gs.gata2 <- read.delim(file.gata2,header=F)
pred.gata2 <- pred.gs.gata2[,2]
gs.gata2 <- pred.gs.gata2[,3]
names(gs.gata2) <- pred.gs.gata2[,1]
names(pred.gata2) <- names(gs.gata2)
#pred.gata2 <- sort(pred.gata2,decreasing=T)
gs.gata2 <- gs.gata2>0
AUPR.res.gata2 <- calcAupr(pred.gata2,gs.gata2)

file.yy1 <- "~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_cent_format"
pred.gs.yy1 <- read.delim(file.yy1,header=F)
pred.yy1 <- pred.gs.yy1[,3]
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
AUPR.res.yy1.cent <- calcAupr(pred.yy1,gs.yy1)

file.gata2 <- "~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_cent_format"
pred.gs.gata2 <- read.delim(file.gata2,header=F)
pred.gata2 <- pred.gs.gata2[,3]
gs.gata2 <- pred.gs.gata2[,4]
names(gs.gata2) <- paste(pred.gs.gata2[,1],pred.gs.gata2[,2],sep= ":")
names(pred.gata2) <- names(gs.gata2)
gs.gata2 <- gs.gata2>0
AUPR.res.gata2.cent <- calcAupr(pred.gata2,gs.gata2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functional enrichment

file.yy1 <- "~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_htseq_format_cut"
pred.gs.yy1 <- read.delim(file.yy1,header=F)
pred.yy1 <- pred.gs.yy1[,3]
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
pred.yy1 <- sort(pred.yy1,decreasing=T)
gene.list.yy1 <- names(pred.yy1)
correl.vec.yy1 <- pred.yy1
gene.set.yy1 <- names(gs.yy1)[which(gs.yy1 > 0)]

GSEA.res.yy1 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.yy1, gene.set=gene.set.yy1,weighted.score.type = 1, correl.vector = correl.vec.yy1)

file.gata2 <- "~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_htseq_format_cut"
pred.gs.gata2 <- read.delim(file.gata2,header=F)
pred.gata2 <- pred.gs.gata2[,3]
gs.gata2 <- pred.gs.gata2[,4]
names(gs.gata2) <- paste(pred.gs.gata2[,1],pred.gs.gata2[,2],sep= ":")
names(pred.gata2) <- names(gs.gata2)
pred.gata2 <- sort(pred.gata2,decreasing=T)
gene.list.gata2 <- names(pred.gata2)
correl.vec.gata2 <- pred.gata2
gene.set.gata2 <- names(gs.gata2)[which(gs.gata2 > 0)]

GSEA.res.gata2 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.gata2, gene.set=gene.set.gata2,weighted.score.type = 1, correl.vector = correl.vec.gata2)

#file.max <- "~/Documents/google-python-exercises/gold_standard/luciferase_MAX_hg19_htseq_format_cut"
#pred.gs.max <- read.delim(file.max,header=F)
#pred.max <- pred.gs.max[,3]
#gs.max <- pred.gs.max[,4]
#names(gs.max) <- paste(pred.gs.max[,1],pred.gs.max[,2],sep= ":")
#names(pred.max) <- names(gs.max)
#pred.max <- sort(pred.max,decreasing=T)
#gene.list.max <- names(pred.max)
#correl.vec.max <- pred.max
#gene.set.max <- names(gs.max)[which(gs.max > 0)]

#GSEA.res.max <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.max, gene.set=gene.set.max,weighted.score.type = 1, correl.vector = correl.vec.max)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functional enrichment ks test
file.yy1 <- "~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_cent_format"
pred.gs.yy1 <- read.delim(file.yy1,header=F)
pred.yy1 <- pred.gs.yy1[,3]
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
pred.yy1 <- sort(pred.yy1,decreasing=T)
gene.list.yy1 <- names(pred.yy1)
correl.vec.yy1 <- pred.yy1
gene.set.yy1 <- names(gs.yy1)[which(gs.yy1 > 0)]

GSEA.res.yy1.cent <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.yy1, gene.set=gene.set.yy1,weighted.score.type = 0)
 correl.vector = correl.vec.yy1

file.gata2 <- "~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_cent_format"
pred.gs.gata2 <- read.delim(file.gata2,header=F)
pred.gata2 <- pred.gs.gata2[,3]
gs.gata2 <- pred.gs.gata2[,4]
names(gs.gata2) <- paste(pred.gs.gata2[,1],pred.gs.gata2[,2],sep= ":")
names(pred.gata2) <- names(gs.gata2)
pred.gata2 <- sort(pred.gata2,decreasing=T)
gene.list.gata2 <- names(pred.gata2)
correl.vec.gata2 <- pred.gata2
gene.set.gata2 <- names(gs.gata2)[which(gs.gata2 > 0)]

GSEA.res.gata2.cent <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.gata2, gene.set=gene.set.gata2,weighted.score.type = 0)
#correl.vector = correl.vec.gata2

file.max <- "~/data/gold_standard/Centipede_K562_hg19_MAX_format"
pred.gs.max <- read.delim(file.max,header=F)
pred.max <- pred.gs.max[,3]
gs.max <- pred.gs.max[,4]
names(gs.max) <- paste(pred.gs.max[,1],pred.gs.max[,2],sep= ":")
names(pred.max) <- names(gs.max)
pred.max <- sort(pred.max,decreasing=T)
gene.list.max <- names(pred.max)
correl.vec.max <- pred.max
gene.set.max <- names(gs.max)[which(gs.max > 0)]

GSEA.res.max <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.max, gene.set=gene.set.max,weighted.score.type = 0, correl.vector = correl.vec.max)

within 50bp of centipede start and end
GATA2100_DGF_K562_htseq_output_new.txt ~/data/gold_standard/Centipede_K562_hg19_GATA2.bed
269
565
877047
xc406@steve:~/code$ python htseqCent.py ~/data/htseq_output/MAX100_DGF_K562_htseq_output_new.txt ~/data/gold_standard/Centipede_K562_hg19_MAX.bed
384
708
1207814
xc406@steve:~/code$ python htseqCent.py ~/data/htseq_output/YY1100_DGF_K562_htseq_output_new.txt ~/data/gold_standard/Centipede_K562_hg19_YY1.bed
1091
2180
1587311

within 10bp of centipede start
gata2
number of gs 269
TP 399
number of total predictions 877047
max
number of gs 384
TP 428
number of total predictions 1207814
yy1
number of gs 1091
TP 1378
number of total predictions 1587311

within 10bp of centipede stand and half-site
gata2
number of gs 269
TP 393
number of total predictions 877047
max
number of gs 384
TP 412
number of total predictions 1207814
yy1
number of gs 1091
TP 1256
number of total predictions 1587311

file.yy1 <- "~/Documents/google-python-exercises/gold_standard/Centipede_K562_hg19_YY1_format"
pred.gs.yy1 <- read.delim(file.yy1,header=F)
pred.yy1 <- pred.gs.yy1[,3]
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
pred.yy1 <- sort(pred.yy1,decreasing=T)
gene.list.yy1 <- names(pred.yy1)
correl.vec.yy1 <- pred.yy1
gene.set.yy1 <- names(gs.yy1)[which(gs.yy1 > 0)]

GSEA.res.yy1 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.yy1, gene.set=gene.set.yy1,weighted.score.type = 1, correl.vector = correl.vec.yy1)

file.gata2 <- "~/Documents/google-python-exercises/gold_standard/Centipede_K562_hg19_GATA2_format"
pred.gs.gata2 <- read.delim(file.gata2,header=F)
pred.gata2 <- pred.gs.gata2[,3]
gs.gata2 <- pred.gs.gata2[,4]
names(gs.gata2) <- paste(pred.gs.gata2[,1],pred.gs.gata2[,2],sep= ":")
names(pred.gata2) <- names(gs.gata2)
pred.gata2 <- sort(pred.gata2,decreasing=T)
gene.list.gata2 <- names(pred.gata2)
correl.vec.gata2 <- pred.gata2
gene.set.gata2 <- names(gs.gata2)[which(gs.gata2 > 0)]

GSEA.res.gata2 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.gata2, gene.set=gene.set.gata2,weighted.score.type = 1, correl.vector = correl.vec.gata2)

file.max <- "~/Documents/google-python-exercises/gold_standard/Centipede_K562_hg19_MAX_format"
pred.gs.max <- read.delim(file.max,header=F)
pred.max <- pred.gs.max[,3]
gs.max <- pred.gs.max[,4]
names(gs.max) <- paste(pred.gs.max[,1],pred.gs.max[,2],sep= ":")
names(pred.max) <- names(gs.max)
pred.max <- sort(pred.max,decreasing=T)
gene.list.max <- names(pred.max)
correl.vec.max <- pred.max
gene.set.max <- names(gs.max)[which(gs.max > 0)]

GSEA.res.max <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.max, gene.set=gene.set.max,weighted.score.type = 1, correl.vector = correl.vec.max)

gns <- length(GSEA.res.yy1$RES)
GSEA.res.perm.ss <- numeric(length(gns))
GSEA.res.perm.s <- numeric(length(gns))
GSEA.ES.res.perm <- numeric()
n=0
for (i in 1:N.perm){
	x <- GSEA.EnrichmentScore(run.perm=T,gene.list=gene.list.yy1, gene.set=gene.set.yy1,weighted.score.type = 1, correl.vector = correl.vec.yy1)
	if(x$ES > 0){
		n=n+1
		GSEA.ES.res.perm <- c(GSEA.ES.res.perm,x$ES)
		res.perm <- x$RES
		GSEA.res.perm.ss <- GSEA.res.perm.ss + res.perm^2
		GSEA.res.perm.s <- GSEA.res.perm.s + res.perm
	}
}

GSEA.res.perm.std <- sqrt( (GSEA.res.perm.ss - ((GSEA.res.perm.s^2)/n)) / (n-1) )
GSEA.res.perm.avg <- GSEA.res.perm.s / n
GSEA.res.yy1$pval <- max(length(which(GSEA.ES.res.perm >= GSEA.res.yy1$ES))/N.perm,1/N.perm)
#GSEA.res.bbz$pval <- max(length(which(GSEA.ES.res.perm >= GSEA.res.bbz$ES))/N.perm,1/N.perm)
#GSEA.res.icos$pval <- max(length(which(GSEA.ES.res.perm >= GSEA.res.icos$ES))/N.perm,1/N.perm)


rand.perm.yy1 <- GSEA.EnrichmentScore(run.perm=T,gene.list=gene.list.yy1, gene.set=gene.set.yy1,weighted.score.type = 1, correl.vector = correl.vec.yy1)
rand.perm.gata2 <- GSEA.EnrichmentScore(run.perm=T,gene.list=gene.list.gata2, gene.set=gene.set.gata2,weighted.score.type = 1, correl.vector = correl.vec.gata2)
rand.perm.max <- GSEA.EnrichmentScore(run.perm=T,gene.list=gene.list.max, gene.set=gene.set.max,weighted.score.type = 1, correl.vector = correl.vec.max)
res <- list(GSEA.res.yy1,GSEA.res.gata2,GSEA.res.max,rand.perm.yy1,rand.perm.gata2,rand.perm.max)
#res <- list(GSEA.res.yy1,GSEA.res.gata2,GSEA.res.yy1.cent,GSEA.res.gata2.cent)
names(res) <- c("YY1","GATA2","MAX","Random")
cls <- c("firebrick3",rgb(255, 106, 106,maxColorValue=255, 150),rgb(154, 205, 50,maxColorValue=255, 150),"lavender")


coords <- list()
coords[[6]] <- c(-1,-0.8)
coords[[5]] <- c(-0.7,-0.5)
coords[[4]] <- c(-0.4,-0.2)
coords[[3]] <- c(-0.1,0.1)
coords[[2]] <- c(0.2,0.4)
coords[[1]] <- c(0.5,0.7)
n <- 10
for(k in 1:length(res)){
    l <- length(res[[k]]$RES)
	if(k ==1){
		plot(y=c(-0.05,0.05),x=c(1,1),col="whitesmoke",type="l",lty=1,ylim=c(0,n),lwd=0.1, xlab="",ylab="",axes=F,xlim=c(-1,1))
	}	
	y <- c(0,0,n,n)
	x <- c(coords[[k]],rev(coords[[k]]))
	polygon(x,y, col="whitesmoke")
	s <- res[[k]]$indicator
	ix <- which(s == 1)
	for(i in ix){
		lines(x=coords[[k]],y=c(n-i*n/l,n-i*n/l),col="firebrick1",lwd=.1)
	}
}
for(k in 1:length(res)){
    l <- length(res[[k]]$RES)
	i = res[[k]]$arg.ES
	x <- c(coords[[k]])
	lines(x=x,y=c(n-i*n/l,n-i*n/l),col="royalblue",lwd=3)
}

plot(x=1:length(res[[1]]$RES),y=res[[1]]$RES,type="l",lwd=3,col='salmon',cex.lab=1,ylim=c(-0.2,1),ylab="Enrichment Score (ES)",xaxt="n",xlab="",cex.axis=1,las=1)
lines(x=1:length(res[[2]]$RES),y=res[[2]]$RES,lwd=3,col='rosybrown1',cex.axis=1,las=1)
lines(x=1:length(rand.perm.yy1$RES),y=rand.perm.yy1$RES,lwd=3,col='royalblue',cex.axis=1,las=1)
lines(x=1:length(rand.perm.gata2$RES),y=rand.perm.gata2$RES,lwd=3,col='lightblue',cex.axis=1,las=1)
legend("topright",lty=1,lwd=2,legend=c(paste(sep="","ES.YY1=",0.66),paste(sep="","ES.GATA2=", 0.77), paste(sep="","ES.YY1.Random=",0.35 ),
paste(sep="","ES.GATA2.Random=", 0.42)),col=c("salmon","rosybrown","royalblue","lightblue"))

########################################################################################################################################################################################

color = c('skyblue','skyblue','skyblue','skyblue','blue','blue','blue','blue','blue','blue','darkblue','darkblue','darkblue','darkblue','darkblue')
i <- 1
plot(x= c(i-1,i),y=c(aupr.avg[i],aupr.avg[i]),type = 'l',xaxt = 'n',col = color[i],xlim = c(0,15),ylim=c(0,1),lwd=1.5,ylab = ('Area Under the Precision Recall Curve'), xlab = (''))
for (i in 2:5){
        lines(x= c(i-1,i),y=c(aupr.avg[i],aupr.avg[i]),xaxt = 'n',type = 'l',col = color[i],xlim = c(0,5),ylim=c(0,1),lwd=1.5)
}
aupr.avg <- (aupr.rorc[1,]+aupr.stat3[1,]+aupr.irf4[1,])/3
barplot(aupr.avg,col = c('skyblue','skyblue','skyblue','skyblue','blue','blue','blue','blue','blue','blue','darkblue','darkblue','darkblue','darkblue','darkblue'),xlim = c(0,15),ylim=c(0,1),space=0,names.arg=c("","","","","","","","","","","","","","",""),ylab = ('Area Under the Precision Recall Curve'),width=1)

points(1:5 - 0.5, aupr.irf4[1,],col = 'red',pch=22,xlim =c(0,5),ylim = c(0,1),type = 'p',lwd=2)
points(1:5 - 0.5, aupr.rorc[1,],col = 'violet',pch=24,xlim =c(0,5),ylim = c(0,1),type = 'p',lwd=2)
points(1:5 - 0.5, aupr.stat3[1,],col = 'orange',pch=23,xlim =c(0,5),ylim = c(0,1),type = 'p',lwd=2)

legend("topleft",pch = c(24,23,22),col=c('violet','orange','red'),legend=c('RORC','STAT3','IRF4'),bty='n',pt.lwd=2)