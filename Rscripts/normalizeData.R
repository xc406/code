##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Aug 2011 nanomed project (Mike Milones human cd4/cd8 data)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
x <- unlist(strsplit(date()," +",perl=TRUE))
date.is <- paste(x[2],x[3],x[5],sep="_")

setwd("/Users/xichen/")
library(affy)

#set paths
path.input <- "data/Microarray_1/cel_files/"
path.output <- "data/Microarray_1/output/"
# get cel files
cel.files <- list.files(path.input)

# normalize data (probeset level)
eset <- justRMA(filenames=paste(sep="",path.input,cel.files))
#save(eset,file=paste(sep="",path.output,"nanomed_and_dustin_data_este.RData"))
d <- exprs(eset)
#save(d,file=paste(sep="",path.output,"nanomed_and_dustin_data_matrix.RData"))

# convert probeset to gene names
library("hugene11sttranscriptcluster.db")
gMap <- hugene11sttranscriptclusterSYMBOL
mProbes <- mappedkeys(gMap)
probNames <- rownames(d)
gMap1 <- as.list(gMap[mProbes]) #list that contains affy_pro
#I assume the mappedProbes are sorted (by index) 
#  -from what I can tell, this assumption is valid
mappedProbes <- sort(which( (probNames %in% names(gMap1)) == T))
any( mappedProbes != sort(mappedProbes)) # FALSE
notMappedProbes <- which( (probNames %in% names(gMap1)) == F)
mappedNames <- as.character(gMap1[ which(names(gMap1) %in% rownames(d)[mappedProbes]) ])
rownames(d)[mappedProbes] <- mappedNames
d <- d[-notMappedProbes,,drop=F] 

# get unique genes matrix by taking mean expression of gene names that appear more than one time
gn.nms <- rownames(d) <- toupper(rownames(d))
d.unique <- matrix(0, nrow = length(unique(gn.nms)), ncol = ncol(d))
seen.gns <- character()
cnt <- 0
for (i in 1:dim(d)[1]){
	cr.gn <- gn.nms[i]
	if (cr.gn %in% seen.gns){
		next	
	} else{
		seen.gns <- c(seen.gns,cr.gn)
		ix <- which(gn.nms %in% cr.gn)
		if(length(ix)==0) {
			stop("length(ix) = 0 at i=",ix[1]," gene name is:", cr.gn)
		}
		if(length(ix)>1){
			cnt <- cnt + 1
			d.unique[cnt,] <- colMeans(d[ix,])
		} else {
			cnt <- cnt + 1
			d.unique[cnt,] <- d[ix,]
		}
	}
}

rownames(d.unique) <- seen.gns
# make pretty colnames
 x <- sapply(strsplit(colnames(d),'_',fixed = T), function(i) i[2])
 w <- read.table(sep="\t",header=T,as.is=T,file=paste(sep="",path.input,"../mapping.txt"))
 x.nm.map <- character(length(x))
 for(i in 1:length(x)){
   ix <- which(w[,1]==x[i])
   x.nm.map[ix] <- w[ix,2]
 }
 colnames(d.unique) <- x.nm.map

d <- d.unique
write.table(d.unique,file=paste(path.output,"sonia_28_icos_bb_time_series_data_matrix","_",date.is,".xls",sep=""),sep="\t")
save(d,file=paste(sep="",path.output,"sonia_28_icos_bb_time_series_data_matrix","_",date.is,".RData"))
# produce quality ctrl dendrogram
hc <- hclust(dist(t(d)), "ave")
boxplot(d)

# postscript looks better with the small fonts
pdf(paste(sep="",path.output,"experiments_dendrogram_",date.is,".pdf"),fonts="Times")
	plot(hc,cex=.75,lwd=.5)
dev.off()

##grep("inf",rownames(d),ignore.case=T,value=T)












