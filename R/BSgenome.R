library(BSgenome.Mmusculus.UCSC.mm9)

cluster <- read.delim("~/data/cluster_background10.txt", header = F)
cluster.df <- data.frame(chr = cluster[,1], start = cluster[,2]-100, end = cluster[,3]+100, id = cluster[,4])
#names(cluster.df) <- c("chr", "start", "end", "id")
x <- getSeq(Mmusculus, cluster.df$chr,cluster.df$start,cluster.df$end)
temp_file <- file.path("~/data/", "temp.fa")
n <- as.character(cluster.df$id)
writeFASTA(x, file=temp_file, desc = n)

