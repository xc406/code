library(BSgenome.Hsapiens.UCSC.hg19)

##write genomic regions in bed format to fasta
cluster <- read.delim("~/data/hg19genesplus10kb.bed", header = F)
cluster.df <- data.frame(chr = cluster[,1], start = cluster[,2], end = cluster[,3], gname = cluster[,4])
#names(cluster.df) <- c("chr", "start", "end", "id")
x <- getSeq(hsapiens, cluster.df$chr,cluster.df$start,cluster.df$end)
temp_file <- file.path("~/data/", "temp.fa")
n <- as.character(cluster.df$gname)
writeFASTA(x, file=temp_file, desc = n)