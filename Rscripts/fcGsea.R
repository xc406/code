rm(list=ls())
setwd("/Users/xichen/")

# source required functions
source("code/microarray_r_scripts/util.R")
library(samr)

# set paths
path.input <- "data/Microarray_1/output/"
path.output <- "data/Microarray_1/results/"

# set file names
file.data <- paste(sep="",path.input,"sonia_28_icos_bb_time_series_data_matrix_Dec_11_2012.RData")
file.tfs <- paste(sep="",path.input,"humanTFs.RData")

# load data (variable name is d)
load(file.data)
# load human tf names humanTFNames
load(file.tfs)

# load literature curated Th17 list
# 74 genes as literature gs in mm9
gs <- as.character(read.delim("~/data/Microarray_1/output/gshg.txt", header=F)$V1)
tmp<-c(rep(1,length(gs)))
names(tmp) <- toupper(gs)
gs <- tmp

gene.set <- names(gs)[which(gs == TRUE)]
source('~/code/microarray_r_scripts/GSEA.R')

# foldchange for three CARs per donor
# 4hr vs 0hr
pred.icos.a <- 2^d[,12]/2^d[,3]
pred.icos.b <- 2^d[,15]/2^d[,6]
pred.icos.c <- 2^d[,18]/2^d[,9]
# 4days vs 0hr
#pred.icos.a <- 2^d[,39]/2^d[,3]
#pred.icos.b <- 2^d[,42]/2^d[,6]
#pred.icos.c <- 2^d[,45]/2^d[,9]

pred.icos.mean <- (pred.icos.a + pred.icos.b + pred.icos.c)/3
pred.icos.mean.sort <- sort(pred.icos.mean,decreasing=T)
gene.list.icos.mean <- names(pred.icos.mean.sort) ## order gene names based on mean foldchange over three donors for GSEA test

pred.bbz.a <- 2^d[,11]/2^d[,2]
pred.bbz.b <- 2^d[,14]/2^d[,5]
pred.bbz.c <- 2^d[,17]/2^d[,8]

#pred.bbz.a <- 2^d[,38]/2^d[,2]
#pred.bbz.b <- 2^d[,41]/2^d[,5]
#pred.bbz.c <- 2^d[,44]/2^d[,8]

pred.bbz.mean <- (pred.bbz.a + pred.bbz.b + pred.bbz.c)/3
pred.bbz.mean.sort <- sort(pred.bbz.mean,decreasing=T)
gene.list.bbz.mean <- names(pred.bbz.mean.sort)

pred.z28.a <- 2^d[,10]/2^d[,1]
pred.z28.b <- 2^d[,13]/2^d[,4]
pred.z28.c <- 2^d[,16]/2^d[,7]

#pred.z28.a <- 2^d[,37]/2^d[,1]
#pred.z28.b <- 2^d[,40]/2^d[,4]
#pred.z28.c <- 2^d[,43]/2^d[,7]

pred.z28.mean <- (pred.z28.a + pred.z28.b + pred.z28.c)/3
pred.z28.mean.sort <- sort(pred.z28.mean,decreasing=T)
gene.list.z28.mean <- names(pred.z28.mean.sort)

# sort foldchange values for three CARs per donor
pred.icos.a.sort <- sort(pred.icos.a,decreasing=T)
pred.icos.b.sort <- sort(pred.icos.b,decreasing=T)
pred.icos.c.sort <- sort(pred.icos.c,decreasing=T)

pred.bbz.a.sort <- sort(pred.bbz.a,decreasing=T)
pred.bbz.b.sort <- sort(pred.bbz.b,decreasing=T)
pred.bbz.c.sort <- sort(pred.bbz.c,decreasing=T)

pred.z28.a.sort <- sort(pred.z28.a,decreasing=T)
pred.z28.b.sort <- sort(pred.z28.b,decreasing=T)
pred.z28.c.sort <- sort(pred.z28.c,decreasing=T)

#get vecs for gsea
gene.list.icos.a <- names(pred.icos.a.sort) # names
gene.list.icos.b <- names(pred.icos.b.sort)
gene.list.icos.c <- names(pred.icos.c.sort)

gene.list.bbz.a <- names(pred.bbz.a.sort)
gene.list.bbz.b <- names(pred.bbz.b.sort)
gene.list.bbz.c <- names(pred.bbz.c.sort)

gene.list.z28.a <- names(pred.z28.a.sort)
gene.list.z28.b <- names(pred.z28.b.sort)
gene.list.z28.c <- names(pred.z28.c.sort)

correl.vec.icos.a <- pred.icos.a.sort # values
correl.vec.icos.b <- pred.icos.b.sort
correl.vec.icos.c <- pred.icos.c.sort

correl.vec.bbz.a <- pred.bbz.a.sort
correl.vec.bbz.b <- pred.bbz.b.sort
correl.vec.bbz.c <- pred.bbz.c.sort

correl.vec.z28.a <- pred.z28.a.sort
correl.vec.z28.b <- pred.z28.b.sort
correl.vec.z28.c <- pred.z28.c.sort

# get a list of all genes
gns <- names(pred.icos.a.sort)

GSEA.res.icos <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.icos.mean, gene.set=gene.set,weighted.score.type = 1, correl.vector = pred.icos.mean.sort)
GSEA.res.bbz <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.bbz.mean, gene.set=gene.set,weighted.score.type = 1, correl.vector = pred.bbz.mean.sort)
GSEA.res.z28 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.z28.mean, gene.set=gene.set,weighted.score.type = 1, correl.vector = pred.z28.mean.sort)

GSEA.res.icos.a <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.icos.a, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos.a)
GSEA.res.icos.b <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.icos.b, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos.b)
GSEA.res.icos.c <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.icos.c, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos.c)

GSEA.res.bbz.a <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.bbz.a, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.bbz.a)
GSEA.res.bbz.b <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.bbz.b, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.bbz.b)
GSEA.res.bbz.c <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.bbz.c, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.bbz.c)

GSEA.res.z28.a <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.z28.a, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.z28.a)
GSEA.res.z28.b <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.z28.b, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.z28.b)
GSEA.res.z28.c <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.z28.c, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.z28.c)

# get permutations for gsea analysis
GSEA.res.perm.ss <- numeric(length(gns))
GSEA.res.perm.s <- numeric(length(gns))
GSEA.ES.res.perm <- numeric()

n=0
N.perm <- 10000 
for (i in 1:N.perm){
	x <- GSEA.EnrichmentScore(run.perm=T,gene.list=gene.list.icos.a, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos.a)
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

# sort GSEA data by gene names to calculate the mean and std
#GSEA.res.icos.a.sort <- GSEA.res.icos.a$RES[sort(names(GSEA.res.icos.a$RES))]
#GSEA.res.icos.b.sort <- GSEA.res.icos.b$RES[sort(names(GSEA.res.icos.b$RES))]
#GSEA.res.icos.c.sort <- GSEA.res.icos.c$RES[sort(names(GSEA.res.icos.c$RES))]

#GSEA.res.bbz.a.sort <- GSEA.res.bbz.a$RES[sort(names(GSEA.res.bbz.a$RES))]
#GSEA.res.bbz.b.sort <- GSEA.res.bbz.b$RES[sort(names(GSEA.res.bbz.b$RES))]
#GSEA.res.bbz.c.sort <- GSEA.res.bbz.c$RES[sort(names(GSEA.res.bbz.c$RES))]

#GSEA.res.z28.a.sort <- GSEA.res.z28.a$RES[sort(names(GSEA.res.z28.a$RES))]
#GSEA.res.z28.b.sort <- GSEA.res.z28.b$RES[sort(names(GSEA.res.z28.b$RES))]
#GSEA.res.z28.c.sort <- GSEA.res.z28.c$RES[sort(names(GSEA.res.z28.c$RES))]

#GSEA.res.icos.mean <- (GSEA.res.icos.a.sort+GSEA.res.icos.b.sort+GSEA.res.icos.c.sort)/3
#GSEA.res.bbz.mean <- (GSEA.res.bbz.a.sort+GSEA.res.bbz.b.sort+GSEA.res.bbz.c.sort)/3
#GSEA.res.z28.mean <- (GSEA.res.z28.a.sort+GSEA.res.z28.b.sort+GSEA.res.z28.c.sort)/3

GSEA.res.icos.mean <- (as.matrix(GSEA.res.icos.a$RES)+as.matrix(GSEA.res.icos.b$RES)+as.matrix(GSEA.res.icos.c$RES))/3
GSEA.res.bbz.mean <- (as.matrix(GSEA.res.bbz.a$RES)+as.matrix(GSEA.res.bbz.b$RES)+as.matrix(GSEA.res.bbz.c$RES))/3
GSEA.res.z28.mean <- (as.matrix(GSEA.res.z28.a$RES)+as.matrix(GSEA.res.z28.b$RES)+as.matrix(GSEA.res.z28.c$RES))/3

GSEA.res.icos.ss <- as.matrix(GSEA.res.icos.a$RES^2)+as.matrix(GSEA.res.icos.b$RES^2)+as.matrix(GSEA.res.icos.c$RES^2)
GSEA.res.bbz.ss <- as.matrix(GSEA.res.bbz.a$RES^2)+as.matrix(GSEA.res.bbz.b$RES^2)+as.matrix(GSEA.res.bbz.c$RES^2)
GSEA.res.z28.ss <- as.matrix(GSEA.res.z28.a$RES^2)+as.matrix(GSEA.res.z28.b$RES^2)+as.matrix(GSEA.res.z28.c$RES^2)

GSEA.res.icos.std <-   sqrt(( GSEA.res.icos.ss - (((GSEA.res.icos.mean *3) ^2) /3) )  / (3-1) )
GSEA.res.bbz.std <- sqrt((GSEA.res.bbz.ss - (((GSEA.res.bbz.mean *3) ^2) /3)) / (3-1) )
GSEA.res.z28.std <- sqrt((GSEA.res.z28.ss - (((GSEA.res.z28.mean *3) ^2) /3)) / (3-1) )

#GSEA.res.icos.mean <- GSEA.res.icos.mean[names(GSEA.res.icos$RES)]
#GSEA.res.bbz.mean <- GSEA.res.bbz.mean[names(GSEA.res.bbz$RES)]
#GSEA.res.z28.mean <- GSEA.res.z28.mean[names(GSEA.res.z28$RES)]

#GSEA.res.icos.std <- GSEA.res.icos.std[names(GSEA.res.icos$RES)]
#GSEA.res.bbz.std <- GSEA.res.bbz.std[names(GSEA.res.bbz$RES)]
#GSEA.res.z28.std <- GSEA.res.z28.std[names(GSEA.res.z28$RES)]

res <- list(GSEA.res.icos.mean,(GSEA.res.icos.mean+GSEA.res.icos.std),(GSEA.res.icos.mean-GSEA.res.icos.std),GSEA.res.bbz.mean,(GSEA.res.bbz.mean+GSEA.res.bbz.std),(GSEA.res.bbz.mean-GSEA.res.bbz.std),
	GSEA.res.z28.mean,(GSEA.res.z28.mean+GSEA.res.z28.std),(GSEA.res.z28.mean-GSEA.res.z28.std),GSEA.res.perm.avg,(GSEA.res.perm.std+GSEA.res.perm.avg),(GSEA.res.perm.avg-GSEA.res.perm.std))

cls <- c("red",rgb(255, 106, 106,maxColorValue=255, 150),rgb(255, 106, 106,maxColorValue=255, 150),"black", "darkgray","darkgray","grey","lightgray","lightgray","blue","lightcyan","lightcyan")
nf <- layout(rbind(1,2),heights=c(3,1))
op <- par(mar=c(0,7,0.5,0.5),mgp=c(4, 1, 0),font=2,font.axis=2,font.lab=2)

# ES plot

	n <- length(res[[1]])  
	x <- 1:n	
	i <- 1
		plot(x=x,y=res[[i]],lwd=3,col=cls[i],cex.axis=2,las=1,type="l",cex.lab=2,ylim=c(-0.2,1),ylab="Enrichment Score (ES)",xaxt="n",xlab="",cex.axis=1.5,las=1)
		#lines(x=x,y=res[[i+1]],lwd=3,col=cls[i+1],cex.axis=2,las=1)
		#lines(x=x,y=res[[i+2]],lwd=3,col=cls[i+2],cex.axis=2,las=1)
		polygon(x=c(x),y=c(res[[i+1]]), col=cls[i+1],border=NA)
		polygon(x=c(x),y=c(res[[i+2]]), col="white",border=NA)
	i <- 7	
		polygon(x=c(x),y=c(res[[i+1]]), col=cls[i+1],border=NA)
		lines(x=x,y=res[[i]],lwd=3,col=cls[i],cex.axis=2,las=1)
		polygon(x=c(x),y=c(res[[i+2]]), col="white",border=NA)		
	i <- 4	
		polygon(x=c(x),y=c(res[[i+1]]), col=cls[i+1],border=NA)
		lines(x=x,y=res[[i]],lwd=3,col=cls[i],cex.axis=2,las=1)		
		polygon(x=c(x),y=c(res[[i+2]]), col="white",border=NA)

		#polygon(x=c(x),y=c(res[[i]]), col=cls[i+1],border=NA)
		#polygon(x=c(x),y=c(res[[i+2]]), col="white",border=NA)	
		
	i <- 10	
		polygon(x=c(x),y=c(res[[i+1]]), col=cls[i+1],border=NA)		
		polygon(x=c(x,rev(x)),y=c(res[[i+1]],rev(res[[i+2]])),col=cls[i+1],border=NA)
		lines(x=x,y=res[[i]],type="l",lwd=3,col=cls[i])	
		
legend(y=1,x=(n-6000),legend=c("28z","BBz","ICOSz","Random"),lty=1,lwd=3,col=c(cls[7],cls[4],cls[1],cls[10]),cex=1.2,bty = "n")

res <- list(GSEA.res.icos,GSEA.res.bbz,GSEA.res.z28)
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
		lines(y=coords[[k]],x=c(i,i),col="darkred",lwd=1)
	}
}
for(k in 1:length(res)){
	i = res[[k]]$arg.ES
	y <- c(coords[[k]])
	lines(y=y,x=c(i,i),col="darkblue",lwd=3)
}

mtext("ICOSz",side=2,cex=1.5,las=1,line = -0.5,at=0.75)
mtext("BBz",side=2,cex=1.5,las=1,line = -0.5,at=0)
mtext("28z",side=2,cex=1.5,las=1,line = -0.5,at=-0.75)
