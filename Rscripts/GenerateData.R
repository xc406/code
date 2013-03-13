rm(list=ls())
## gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"
load("/Users/xichen/data/results_combine_data_exprsn_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1.Rdata")
#pred.file <- 
#knockout_rna_zscores.all <- read.delim(gs.file)
k_r_network_stat3 <- res[["activation"]][["STAT3"]][["th17"]][,"KRI"] ##knockout_rna_zscores.all[,"IRF4"]
#names(k_r_network_stat3) <- toupper(rownames(knockout_rna_zscores.all))
k_r_network_stat3 <- k_r_network_stat3[which(k_r_network_stat3>0)]

#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/htseq_output/final/Irf4_htseq_output_avg", header=F)
motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/Stat31e3_faire_htseq_output", header=F)
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[-(nrow(motif.and.dnase.reads.stat3.rep1) - 0:4),] #remove meta
#motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>300),2] <- 300
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>0.0), ]
m_d_network_stat3 <- motif.and.dnase.reads.stat3.rep1[,2]
names(m_d_network_stat3) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])
unique.names <- unique(c( names(m_d_network_stat3), names(k_r_network_stat3)))

overlap <- (length(m_d_network_stat3) + length(k_r_network_stat3)) - length(unique.names)
data <- matrix(NA, nr=length(unique.names), nc=2)
rownames(data) <- unique.names
colnames(data) <- c("motif_dnase", "knockout_rna")
for (i in 1:nrow(data)) {
  data[i,1] <- m_d_network_stat3[unique.names[i]]
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
plot(stat3_aupr$rec, stat3_aupr$prec, type = 'lines', main = 'STAT3_KRI', xlab = 'recall', ylab = 'precision')

save.image('~/data/htseqcount.RData')