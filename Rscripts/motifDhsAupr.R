rm(list=ls())
## gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"
load("/Users/xichen/data/results_combine_data_exprsn_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1.Rdata")
#pred.file <- 
#knockout_rna_zscores.all <- read.delim(gs.file)
k_r_network_stat3 <- res[["activation"]][["IRF4"]][["th17"]][,"C"]
#+res[["activation"]][["BATF"]][["th17"]][,"KCRI"]+res[["activation"]][["IRF4"]][["th17"]][,"KCRI"]+res[["activation"]][["MAF"]][["th17"]][,"KCRI"]+res[["activation"]][["RORC"]][["th17"]][,"KCRI"]
k_r_network_stat3 <- k_r_network_stat3[which(k_r_network_stat3>0)]
#k_r_network_stat3 <- read.delim("~/data/Microarray_1/output/gs.txt", header=F)##74 genes as literature gs in mm9

motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/htseq_output/final/IRF4_DGF_htseq_outputmm.txt", header=F)
#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/phylopStat3Window_htseq_output", header=F)
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[-(nrow(motif.and.dnase.reads.stat3.rep1) - 0:4),] #remove meta
#motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>300),2] <- 300
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>2.0), ]


m_d_network_stat3 <- motif.and.dnase.reads.stat3.rep1[,2]
names(m_d_network_stat3) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])
#m_d_network <- m_d_network_stat3[order(m_d_network_stat3, decreasing = T)]
#m_d <- c(rep(NA,200))
#for (i in 1:200){
#	m_d[i] <- m_d_network[i]##motif.and.dnase.reads.stat3.rep1[i,2]
#	}

#m_d <- m_d_network[1:200]

unique.names <- unique(c( names(m_d_network_stat3), names(k_r_network_stat3)))

#overlap <- (length(m_d_network_stat3) + length(k_r_network_stat3)) - length(unique.names)
data <- matrix(NA, nr=length(unique.names), nc=2)
rownames(data) <- unique.names
colnames(data) <- c("motif_dnase", "knockout_rna")
for (i in 1:nrow(data)) {
  #data[i,2] <- m_d[unique.names[i]]
  data[i,1] <- m_d_network_stat3[unique.names[i]]##read count per gene
  data[i,2] <- k_r_network_stat3[unique.names[i]]##kcri network score
  #data[i,2] <- th17_network_stat3[unique.names[i]]##74 th17 genes from literature
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
icosp_aupr <- calcAupr(pred,gs)
plot(icosp_aupr$rec, icosp_aupr$prec, type = 'lines', main = 'IRF4_MOTIF1E3_DGF_VS_KCRI', xlab = 'recall', ylab = 'precision')
#lines(bb_aupr$rec, bb_aupr$prec, col="darkgreen")
#lines(z28_aupr$rec, z28_aupr$prec, col="blue")
legend("topright",lty=1,lwd=2,legend=c("AUPR=0.6955647","AUROC=0.4749878"))
save.image('~/data/htseqcount.RData')
