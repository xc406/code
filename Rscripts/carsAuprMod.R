rm(list=ls())

############ Aviv touchups ##########
use.kcri = T
use.fc = T
num.genes.to.consider.for.gs = 5000
N.perm <- 10000 


# load gold standard genes
if (use.kcri == T){
	# load kcri gold standard genes
	gs <- read.delim("~/data/Microarray_1/output/kcri_corehg.txt", header=F)
	kcri <- gs[,2] ##KCRI network sum score 
	names(kcri) <- toupper(gs[,1])
	kcri <- sort(kcri, decreasing = T)
	kcri.network <- c(rep(NA,num.genes.to.consider.for.gs))
	kcri.network <- kcri[1:num.genes.to.consider.for.gs]
	gs <- kcri.network
	title.tp.2 = 'KCRI network'
} else {
	# load literature curated Th17 list
	##74 genes as literature gs in mm9
	gs <- as.character(read.delim("~/data/Microarray_1/output/gshg.txt", header=F)$V1)
	tmp<-c(rep(1,length(gs)))
	names(tmp) <- toupper(gs)
	gs <- tmp
	title.tp.2 = '74_Th17_genes'
}

# load diff exp resutls
if (use.fc == T){
	diff.exp <- read.delim("~/data/Microarray_1/output/fold_change_4h_vs_ctrl_unpaired.xls", header=T)
	title.tp.1 = "FC"
} else {
	diff.exp <- read.delim("~/data/Microarray_1/output/sam_diff_exp_4h_vs_ctrl_paired.xls", header=T)
	title.tp.1 = "Zscore"
}

diff.exp <- as.matrix(diff.exp)
# for consistancy lets have gene names in upper case form
rownames(diff.exp) <- toupper(rownames(diff.exp))

# now extract icos = "X34h_vs_ctrl"
diff_exp_icos <- diff.exp[,"X34h_vs_ctrl"]##icos 4h vs 0h 
# now extract bbz = "X24h_vs_ctrl"
diff_exp_bbz <- diff.exp[,"X24h_vs_ctrl"]##bbz 4h vs 0h 
# now extract z28 = "X14h_vs_ctrl"
diff_exp_z28 <- diff.exp[,"X14h_vs_ctrl"]##z28 4h vs 0h 

# find genes in gold standard that are not present in microarray data
ix <- which(! names(gs) %in% names(diff_exp_icos))
if(length(ix) > 0){
	gs <- gs[-ix]
}

# get a list of all genes
gns <- rownames(diff.exp)

# create a data matrix to combine all the results with the gold standard
data <- matrix(0, nr=length(gns), nc=4)
rownames(data) <- gns
colnames(data) <- c("gs", "diff_exp_icos","diff_exp_bbz","diff_exp_z28")
data[,"diff_exp_icos"] <- diff_exp_icos 
data[,"diff_exp_bbz"] <- diff_exp_bbz 
data[,"diff_exp_z28"] <- diff_exp_z28 
ix <- which(rownames(data) %in% names(gs))
data[ix,"gs"] <- 1


gs <- data[,1]
gs <- gs>0

pred.icos <- data[,"diff_exp_icos"]
pred.bbz <- data[,"diff_exp_bbz"]
pred.z28 <- data[,"diff_exp_z28"]


pred.icos.sort <- sort(pred.icos,decreasing=T)
pred.bbz.sort <- sort(pred.bbz,decreasing=T)
pred.z28.sort <- sort(pred.z28,decreasing=T)
#get vecs for gsea
gene.list.icos <- names(pred.icos.sort) # names
gene.list.bbz <- names(pred.bbz.sort)
gene.list.z28 <- names(pred.z28.sort)
correl.vec.icos <- pred.icos.sort # values
correl.vec.bbz <- pred.bbz.sort
correl.vec.z28 <- pred.z28.sort

gene.set <- names(gs)[which(gs == TRUE)]
source('~/code/GSEA.R')
GSEA.res.icos <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.icos, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos)
GSEA.res.bbz <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.bbz, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.bbz)
GSEA.res.z28 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.z28, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.z28)
# get permutations for gsea analysis
GSEA.res.perm.ss <- numeric(length(gns))
GSEA.res.perm.s <- numeric(length(gns))
GSEA.ES.res.perm <- numeric()
n=0
for (i in 1:N.perm){
	x <- GSEA.EnrichmentScore(run.perm=T,gene.list=gene.list.icos, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos)
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
GSEA.res.z28$pval <- max(length(which(GSEA.ES.res.perm >= GSEA.res.z28$ES))/N.perm,1/N.perm)
GSEA.res.bbz$pval <- max(length(which(GSEA.ES.res.perm >= GSEA.res.bbz$ES))/N.perm,1/N.perm)
GSEA.res.icos$pval <- max(length(which(GSEA.ES.res.perm >= GSEA.res.icos$ES))/N.perm,1/N.perm)

source('~/code/AUPR/aupr.R')
AUPR.res.icos <- calcAupr(pred.icos,gs)
AUPR.res.bbz <- calcAupr(pred.bbz,gs)
AUPR.res.z28 <- calcAupr(pred.z28,gs)



# AUPR plot
plot.new()
main.title = paste(sep="","CARs_4h_vs_0h_",title.tp.1," vs ", title.tp.2)
plot(AUPR.res.icos$rec, AUPR.res.icos$prec, type = 'l', main = main.title, xlab = 'recall', ylab = 'precision',xlim = c(0,0.3), col="red")
lines(AUPR.res.bbz$rec, AUPR.res.bbz$prec, col="darkgreen")
lines(AUPR.res.z28$rec, AUPR.res.z28$prec, col="blue")
legend("topright",lty=1,lwd=2,legend=c(paste(sep="","AUPR.ICOSz=", round(AUPR.res.icos$AUPR,3)),
	paste(sep="","AUPR.BBz=", round(AUPR.res.bbz$AUPR,3)),paste(sep="","AUPR.28z=", round(AUPR.res.z28$AUPR,3))),col=c("red","darkgreen","blue"))

# AUROC plot
plot.new()
plot(AUPR.res.icos$fpr, AUPR.res.icos$rec, type = 'l', main = main.title, xlab = 'false positive rate', ylab = 'true positive rate', col="red")
lines(AUPR.res.bbz$fpr, AUPR.res.bbz$rec, col="darkgreen")
lines(AUPR.res.z28$fpr, AUPR.res.z28$rec, col="blue")
legend("topright",lty=1,lwd=2,legend=c(paste(sep="","AUROC.ICOSz=", round(AUPR.res.icos$AUROC,3)),paste(sep="","AUROC.BBz=", round(AUPR.res.bbz$AUROC,3)),
	paste(sep="","AUROC.28z=", round(AUPR.res.z28$AUROC,3))),col=c("red","darkgreen","blue"))

res <- list(GSEA.res.icos,GSEA.res.bbz,GSEA.res.z28,GSEA.res.perm.avg,(GSEA.res.perm.std+GSEA.res.perm.avg),(GSEA.res.perm.avg-GSEA.res.perm.std))
names(res) <- c("ICOSz","BBz","28z","Random")
#f.nm <- paste(sep="",path.output,"gsea_",gs.lists[iter],"_",comb.case,".pdf")
#pdf(f.nm)
cls <- c("firebrick3",rgb(255, 106, 106,maxColorValue=255, 150),rgb(154, 205, 50,maxColorValue=255, 150),"lightgray","darkgray", "darkgray")
nf <- layout(rbind(1,2,3),heights=c(3,1,2))
op <- par(mar=c(0,7,0.5,0.5),mgp=c(4, 1, 0))

# ES plot
for(i in 1:length(res)){
	n <- length(res[[1]]$RES)
    #y <- res[[i]]$RES 
	x <- 1:n
	if(i==1){
		plot(x=x,y=res[[i]]$RES,type="l",lwd=3,col=cls[i],cex.lab=2,ylim=c(-0.2,1),ylab="Enrichment Score (ES)",xaxt="n",xlab="",cex.axis=2,las=1)
	} else if (i < 4){
		lines(x=x,y=res[[i]]$RES,lwd=3,col=cls[i],cex.axis=2,las=1)
	} else if (i == 4) {
		lines(x=x,y=res[[i]],lwd=3,col=cls[i],cex.axis=2,las=1)
	} else {
		lines(x=x[seq(from=0,to=length(gns),by=100)],y=res[[i]][seq(from=0,to=length(gns),by=100)],lwd=2,lty=3,col=cls[i],cex.axis=2,las=1)
	}
}
legend(y=1,x=(n-8000),legend=c("ICOSz","BBz","28z","Random"),lty=1,lwd=3,col=cls[1:length(res)],cex=0.8)

coords <- list()
coords[[3]] <- c(-1,-0.5)
coords[[2]] <- c(-0.25,0.25)
coords[[1]] <- c(0.5,1)
for(k in 1:length(res)){
	if(k ==1){
		plot(y=c(-0.05,0.05),x=c(1,1),col="gray89",type="l",lty=1,xlim=c(0,n),lwd=0.1, xlab="",ylab="",axes=F,ylim=c(-1,1))
	}
	x <- c(0,0,n,n)
	y <- c(coords[[k]],rev(coords[[k]]))
	polygon(x,y, col="gray89")
	s <- res[[k]]$indicator
	ix <- which(s == 1)
	for(i in ix){
		lines(y=coords[[k]],x=c(i,i),col="darkred",lwd=2)
	}
}
for(k in 1:length(res)){
	i = res[[k]]$arg.ES
	y <- c(coords[[k]])
	lines(y=y,x=c(i,i),col="darkblue",lwd=5)
}

mtext("ICOSz",side=2,cex=1,las=1,line = -1,at=0.75)
mtext("BBz",side=2,cex=1,las=1,line = -1,at=0)
mtext("28z",side=2,cex=1,las=1,line = -1,at=-0.75)

par(mar=c(4.5,7,0.5,0.5),mgp=c(3, 1, 0))
de <- sort(diff_exp_icos,decreasing=T)
for (i in 1:3){
	if (i==1){
	plot(y=sort(data[,i+1],decreasing=T),x = c(1:length(res[[i]]$RES)),col=cls[i],pch=20,lty=1,xlim=c(0,n),ylim=c(min(de),max(de)),cex.axis=1,cex.lab=1.2,las=1,lwd=3,
		 ylab="FC",xlab="Ranked gene list")

	}else
	lines(x=c(1:length(res[[i]]$RES)),y=sort(data[,i+1],decreasing=T),pch=20,col=cls[i],lwd=3)
}
legend(y=max(de)-1,x=(n-8000),legend=c("ICOSz","BBz","28z"),lty=1,lwd=3,col=cls[c(1,2,3,4)],cex=1)
#lines(x=1:n,y=sort(diff_exp_bbz,decreasing=T),pch=20,col=cls[2],lwd=1)
#lines(x=1:n,y=sort(diff_exp_z28,decreasing=T),pch=20,col=cls[3],lwd=1)
