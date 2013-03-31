######################################################################################
## gsea on motif+dhs --> target prediction
gs <- read.delim("~/data/mycMotifTargets.txt", header=F)##gene set
gs<- gs[,1]


gl <- read.delim("~/data/myc_mcf7_all.txt",header=T)##gene list
gl <- as.matrix(gl)
for (i in 1:length(gl[,3])){
	#if (as.numeric(gl[i,3]) < 0.05){ 
		gl.n[i] <- (-1)*log10(as.numeric(gl[i,3]))#1000
	#}else{
		#gl.n[i] <- (-1)*log10(as.numeric(gl[i,3]))
	#}
}

names(gl.n) <- gl[,1]

ix <- which(! gs %in% names(gl.n))
if(length(ix) > 0){
 	gs <- gs[-ix]
 }

#ix2 <- which(gl.n < 1.30 )
#if(length(ix2) > 0){
#	gl.n <- gl.n[-ix2]
#}

source('~/code/microarray_r_scripts/GSEA.R')
GSEA.res <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=names(gl.n), gene.set=gs,weighted.score.type = 1, correl.vector = gl.n)

# get permutations for gsea analysis
GSEA.res.perm.ss <- numeric(length(gl.n))
GSEA.res.perm.s <- numeric(length(gl.n))
GSEA.ES.res.perm <- numeric()

n=0
N.perm <- 10000 
for (i in 1:N.perm){
	x <- GSEA.EnrichmentScore(run.perm=T,gene.list=names(gl.n), gene.set=gs,weighted.score.type = 1, correl.vector = gl.n)
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

res <- list(GSEA.res$RES,GSEA.res.perm.avg,(GSEA.res.perm.std+GSEA.res.perm.avg),(GSEA.res.perm.avg-GSEA.res.perm.std))
cls <- c("red",rgb(255, 106, 106,maxColorValue=255, 150),"lightgray","lightgray")
nf <- layout(rbind(1,2),heights=c(3,1))
op <- par(mar=c(0,7,0.5,0.5),mgp=c(4, 1, 0))

# ES plot

	n <- length(res[[1]])  
	x <- 1:n
	i <- 1
		plot(x=x,y=res[[i]],type="l",lwd=3,col=cls[i],cex.lab=2,ylim=c(-1,1),ylab="Enrichment Score (ES)",xaxt="n",xlab="",cex.axis=2,las=1)
		

		for (i in 2:4){
			lines(x=x,y=res[[i]],lwd=3,col=cls[i],cex.axis=2,las=1)
		}
		

#######################################################################################