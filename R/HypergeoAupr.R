rm(list=ls())
hypergeo.pvals <- read.delim("~/data/genephyper", header = FALSE)
pvals <- hypergeo.pvals[,2]
names(pvals) <- hypergeo.pvals[,1]
pvals.cutoff <- pvals[which(pvals<0.00001)]
pvals.cutoff.sort <- sort(pvals.cutoff)
# pvals.df <- data.frame(gname = hypergeo.pvals[,1], pvals = hypergeo.pvals[,2])
# pvals.sort.df <- pvals.df[with(pvals.df, order(pvals)), ]
#motif.matrix <- matrix(c(names(pvals), as.numeric(pvals)), nr = length(names(pvals)), nc = 2)
#write.table(motif.matrix, file="~/data/genephyper_sort", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

gwas.file <- "~/Downloads/GWAS-List/gold_standard_celiac_genes_Jun_22_2012.txt"
gs.file <- "~/Downloads/mmc5.xls"
#pred.file <- 
knockout_rna_zscores.all <- read.delim(gs.file)
gwas.cd <- read.delim(gwas.file)
k_r_network_stat3 <- gwas.cd[,"distance"]
k_r_network_stat3 <- knockout_rna_zscores.all[,"TF.sum"]
names(k_r_network_stat3) <- gwas.cd[,"gene_id"]#knockout_rna_zscores.all[,"Gene_id"]#toupper(rownames(knockout_rna_zscores.all))
k_r_network_stat3 <- k_r_network_stat3[which(k_r_network_stat3>0)]

unique.names <- unique(c(names(pvals), names(k_r_network_stat3)))
overlap <- (length(pvals)+length(k_r_network_stat3)) - length(unique.names)

# #motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/htseq_output/final/Irf4_htseq_output_avg", header=F)
# motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/Irf41e3_faire_htseq_output", header=F)
# motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[-(nrow(motif.and.dnase.reads.stat3.rep1) - 0:4),] #remove meta
# #motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>300),2] <- 300
# #motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>0.0), ]
# m_d_network_stat3 <- motif.and.dnase.reads.stat3.rep1[,2]
# names(m_d_network_stat3) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])
# unique.names <- unique(c( names(m_d_network_stat3), names(k_r_network_stat3)))

# overlap <- (length(m_d_network_stat3) + length(k_r_network_stat3)) - length(unique.names)
data <- matrix(NA, nr=length(unique.names), nc=2)
rownames(data) <- unique.names
colnames(data) <- c("motif_dnase", "knockout_rna")
for (i in 1:nrow(data)) {
  data[i,1] <- pvals[unique.names[i]]
  data[i,2] <- k_r_network_stat3[unique.names[i]]
}
#data[which(is.na(data))] <- 0

source('~/code/AUPR/aupr.R')
gs <- data[,2]
gs[which(is.na(gs))] <- 0 

pred <- data[,1]
pred <- pred[which(!is.na(pred))]

# gs <- data[,2]
# colnames(data)
# pred <- data[,1]
# gs>0
gs <- gs>0
# pred <- pred>0
# gs
calcAupr(pred,gs)
stat3_aupr <- calcAupr(pred,gs)
# plot(stat3_aupr$rec, stat3_aupr$prec, type = 'lines', main = 'Irf41e3_faire_htseq_nocutoff', xlab = 'recall', ylab = 'precision')
plot(stat3_aupr$rec, stat3_aupr$prec, type = 'lines', main = 'colocolization_gwas.celiac', xlab = 'recall', ylab = 'precision')
save.image('~/data/htseqcount.RData')