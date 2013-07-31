library(BSgenome.Hsapiens.UCSC.hg19)

##write genomic regions in bed format to fasta
cluster <- read.delim("~/data/hg19ud10.bed", header = F)
cluster.df <- data.frame(chr = cluster[,1], start = cluster[,2], end = cluster[,3], gname = cluster[,4])
#names(cluster.df) <- c("chr", "start", "end", "id")

temp_file <- file.path("~/data/", "temp.fa")
n <- as.character(cluster.df$gname)
x <- getSeq(Hsapiens, cluster.df$chr, cluster.df$start,cluster.df$end)
names(x) <- n
writeXStringSet(x,file = temp_file, format = "fasta")
