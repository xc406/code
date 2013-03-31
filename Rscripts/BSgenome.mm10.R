library(BSgenome.Mmusculus.UCSC.mm10)

##write genomic regions in bed format to fasta
#cluster <- read.delim("~/data/mm9ud5kb.bed", header = F)
#cluster <- read.delim("~/data/AR_coordinates_mm10.txt", header = T)
cluster <- read.delim("~/data/arl_background.txt", header = F)
#cluster.df <- data.frame(chr = cluster[,2], start = as.numeric(cluster[,3]), end = cluster[,4], gname = cluster[,1])
cluster.df <- data.frame(chr = cluster[,1], start = as.numeric(cluster[,2]), end = as.numeric(cluster[,3]), gname = cluster[,4])
#names(cluster.df) <- c("chr", "start", "end", "id")

x <- list()
n <- list()
##split random genomic regions into 10 smaller files
for (i in 1:10){
	x[[i]] <- getSeq(Mmusculus, cluster.df[((i-1)*1000+1):(i*1000),]$chr,cluster.df[((i-1)*1000+1):(i*1000),]$start,cluster.df[((i-1)*1000+1):(i*1000),]$end)
	n[[i]] <- as.character(cluster.df[((i-1)*1000+1):(i*1000),]$gname)
	names(x[[i]]) <- n[[i]]
} 
#x <- getSeq(Mmusculus, cluster.df$chr,cluster.df$start,cluster.df$end)
for (i in 1:10){
	temp_file <- file.path(paste(sep="","~/data/", "c",i,"temp.fa"))
	writeXStringSet(x[[i]],file = temp_file, format = "fasta")
}

#temp_file <- file.path("~/data/", "temp.fa")
#writeXStringSet(x[[1]],file = temp_file, format = "fasta")