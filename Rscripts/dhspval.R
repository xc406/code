## script to assign empirical p-vals to a DHS read counts distribution

x <- read.delim("~/data/MYCud5_DHS1_htseq_output",header = F)
r <- read.delim("~/data/MYCud5_random3_DHS_htseq_output",header = F) ## background read count distribution
ccm <- read.delim("~/data/UniProtID2HGNC_commonname.txt",header = F)
x <- x[-(nrow(x)-0:4),]## strip of the last five rows of meta
r <- r[-(nrow(r)-0:4),]
data <- x[,2]
r.data <- r[,2]
##data <- x2[1:(length(x2)-5)] 
#cdf.func <- ecdf(data)
colnames(x) <- c("gene.name","read.count")
for (i in 1:nrow(x)){
	#x$p.val[i] <- 1-cdf.func(x[i,2]) ## empirical cdf
	x$p.val[i] <- (sum(r.data[which(r.data>x$read.count[i])])+1)/(sum(r.data)+1) ## empirical p-val (s+1)/(n+1)
	
}

#hist(x$read.count,breaks = 1000, xlim = c(0,1500))

## extract genes involved in ccm
ccm.gn <- ccm[,2]

k <- 1
x <- as.matrix(x)
x.ccm <- matrix(NA,nr=length(ccm.gn),nc=3)
for (j in 1:length(ccm.gn)){
	for (i in 1:nrow(x)){
      if (as.character(x[i,1]) == as.character(ccm.gn[j])){
		x.ccm[k,] <- x[i,]
		k <- k+1
		}
	}
}

colnames(x.ccm) <- c("gene.name","read.count","p.val")
x.sort <- x[order(as.numeric(x[,3]),x[,1]),]
x.ccm.sort <- x.ccm[order(as.numeric(x.ccm[,3]),x.ccm[,1]),]
write.table(x.ccm.sort,sep="\t",quote = FALSE,row.names=FALSE,col.names=TRUE,file=paste(sep="",path.input,"myc_mcf7_ccm.txt"))
write.table(x.sort,sep="\t",quote = FALSE,row.names=FALSE,col.names=TRUE,file=paste(sep="",path.input,"myc_mcf7_all.txt"))

length(x.sort[which(as.numeric(x.sort[,3])<0.05)]) ## num of significant targets
length(x.ccm.sort[which(as.numeric(x.ccm.sort[,3])<0.05)]) 
#plot(x = x$read.count, y = x$p.val,xlim=c(0,1000))

