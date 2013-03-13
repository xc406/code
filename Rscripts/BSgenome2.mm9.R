library(BSgenome.Mmusculus.UCSC.mm9)

##write genomic regions in bed format to fasta
cluster <- read.delim("~/data/allgenomicregions", header = F)
cluster.df <- data.frame(chr = cluster[,1], start = cluster[,2], end = cluster[,3], gname = cluster[,4])
#names(cluster.df) <- c("chr", "start", "end", "id")
x <- getSeq(Mmusculus, cluster.df$chr,cluster.df$start,cluster.df$end)
temp_file <- file.path("~/data/", "temp.fa")
n <- as.character(cluster.df$gname)
writeFASTA(x, file=temp_file, desc = n)