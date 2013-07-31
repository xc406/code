library(BSgenome.Mmusculus.UCSC.mm9)

##write genomic regions in bed format to fasta
cluster <- read.delim("~/data/esc/100All_InSilico_dist_id.bed", header = F)
cluster.df <- data.frame(chr = cluster[,1], start = cluster[,2], end = cluster[,3], gname = cluster[,4])
#names(cluster.df) <- c("chr", "start", "end", "id")
x <- getSeq(Mmusculus, cluster.df$chr,cluster.df$start,cluster.df$end)
temp_file <- file.path("~/data/esc/", "esc_insilico_dist.fa")
n <- as.character(cluster.df$gname)
names(x) <- n
writeXStringSet(x,file = temp_file, format = "fasta")
