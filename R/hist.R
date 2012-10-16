#args <- commandArgs()
#gff <- read.table(args[], sep='\t', header=F)
rorc <- read.delim("~/data/htseq_output/final/Stat3_htseq_output_avg", header=F)
rorc <- rorc[-(nrow(rorc) - 0:4),]
#range(rorc[,2])
rorc[which(rorc[,2]>600),2] <- 600
rorc <- rorc[which(rorc[,2]<601.0),]
rorc <- rorc[which(rorc[,2]>0.0),]
#hist(rorc[,2], freq =F, breaks = 500, main = 'Rorc', xlab = 'DHS counts')
hist(rorc[,2], freq =F, breaks = 500, main = 'Stat3', xlab = 'DHS counts')
abline(v=0, col = 2, lty = 1)
