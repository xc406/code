library(MASS)

cluster.count <- read.delim("~/data/mm9cluster1e5fimoout_101112_count", header = F)
background.count <- read.delim("~/data/mm9background1e5fimoout_101212_count", header = F)
cluster.count.df <- data.frame(row.names = cluster.count[,1], cluster.target.count = cluster.count[,2], cluster.tot.count = cluster.count[,3])
background.count.df <- data.frame(row.names = background.count[,1], background.target.count = background.count[,2], background.tot.count = background.count[,3])
# x <- 0:380
all.df <- cbind(cluster.count.df, background.count.df) 
pvals <- rep(NA, length(rownames(all.df)))
names(pvals) <- rownames(all.df)
for (i in rownames(all.df)){ 
  pvals[i] <- phyper(all.df[i,names(all.df)[1]],all.df[i,names(all.df)[3]],38000-all.df[i,names(all.df)[3]],3800, lower.tail = FALSE)
  #names(pvals[i]) <- i
}
length(pvals[(which(pvals<0.00001))])
motif.matrix <- matrix(c(names(pvals), as.numeric(pvals)), nr = length(names(pvals)), nc = 2)
#motif.matrix <- as.matrix(pvals)
#motif.list <- do.call("rbind", as.list(motif.matrix))
motif.matrix.cutoff <- matrix(c(names(pvals[which(pvals<2.2e-16)]), as.numeric(pvals[which(pvals<2.2e-16)])), nr = length(names(pvals[which(pvals<2.2e-16)])), nc = 2)
write.table(motif.matrix, file="~/data/motifphyper.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

pvals.sort <- sort(pvals)
