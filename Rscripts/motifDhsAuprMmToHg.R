rm(list=ls())
## gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"
load("/Users/xichen/data/results_combine_data_exprsn_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1.Rdata")
#pred.file <- 
#knockout_rna_zscores.all <- read.delim(gs.file)

stat3.gs <- res[["activation"]][["STAT3"]][["th17"]][,1:15]
#stat3_gs <- cbind(stat3_kcri,stat3_kc,stat3_kri,stat3_k,stat3_cri,stat3_c,stat3_ri)

mm.to.hg <- read.delim("~/data/HMD_HumanSequence.txt", header=F)
mm.to.hg <- as.matrix(mm.to.hg)


stat3.gs.hg <- matrix(NA,nr = length(intersect(rownames(stat3.gs),toupper(mm.to.hg[,1]))),nc = ncol(stat3.gs))
k <- 1
rownames(stat3.gs.hg) <- rep(NA,length(intersect(rownames(stat3.gs),toupper(mm.to.hg[,1]))))
colnames(stat3.gs.hg) <- colnames(stat3.gs)
for (i in 1:nrow(stat3.gs)){
    for (j in 1:nrow(mm.to.hg)){
	if (toupper(mm.to.hg[j,1]) == rownames(stat3.gs)[i]){
	    stat3.gs.hg[k,] <- stat3.gs[i,]
	    rownames(stat3.gs.hg)[k] <- toupper(mm.to.hg[j,4])
            k <- k+1
	}
    }
}

save.image('~/code/AUPR/th17Stat3MmToHg.RData')

#+res[["activation"]][["BATF"]][["th17"]][,"KCRI"]+res[["activation"]][["IRF4"]][["th17"]][,"KCRI"]+res[["activation"]][["MAF"]][["th17"]][,"KCRI"]+res[["activation"]][["RORC"]][["th17"]][,"KCRI"]

#k_r_network_stat3 <- read.delim("~/data/Microarray_1/output/gs.txt", header=F)##74 genes as literature gs in mm9

motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/htseq_output/final/STAT31e3ud5_DHS_htseq_output", header=F)
#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/phylopStat3Window_htseq_output", header=F)
motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[-(nrow(motif.and.dnase.reads.stat3.rep1) - 0:4),] #remove meta
#motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>300),2] <- 300
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>2.0), ]


m_d_network_stat3 <- motif.and.dnase.reads.stat3.rep1[,2]
names(m_d_network_stat3) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])

unique.names <- unique(c( names(m_d_network_stat3), rownames(stat3.gs.hg)))

#unique.names <- union(names(m_d_network_stat3), rownames(stat3.gs.hg))
#overlap <- (length(m_d_network_stat3) + length(k_r_network_stat3)) - length(unique.names)
data <- matrix(NA, nr=length(unique.names), nc=16)
rownames(data) <- unique.names
colnames(data) <- c("motif_dnase",colnames(stat3.gs.hg))
for (i in 1:nrow(data)) {
  #data[i,2] <- m_d[unique.names[i]]
  data[i,1] <- m_d_network_stat3[unique.names[i]]##read count per gene
  for (j in 2:ncol(data)){
     if (unique.names[i] %in% rownames(stat3.gs.hg)){
        data[i,j] <- stat3.gs.hg[unique.names[i],(j-1)]##kcri network score
     }
  }
  #data[i,2] <- th17_network_stat3[unique.names[i]]##74 th17 genes from literature
}
#data[which(is.na(data))] <- 0

pred <- data[,1]
pred <- pred[which(!is.na(pred))]

source('~/code/AUPR/aupr.R')
aupr_res <- matrix(NA, nc = ncol(stat3.gs.hg), nr = 2)##(2x7)
for (i in 2:ncol(data)){
  gs <- data[,i]
  gs[which(is.na(gs))] <- 0 
  gs <- gs>0
  aupr <- calcAupr(pred,gs)
  aupr_res[,(i-1)] <- c(aupr$AUPR, aupr$AUROC)
}
rownames(aupr_res) <- c("AUPR","AUROC")
colnames(aupr_res) <- colnames(stat3.gs.hg)
aupr_res

#icosp_aupr <- calcAupr(pred,gs)
#plot(icosp_aupr$rec, icosp_aupr$prec, type = 'lines', main = 'IRF4_MOTIF1E3_DGF_VS_KCRI', xlab = 'recall', ylab = 'precision')
#lines(bb_aupr$rec, bb_aupr$prec, col="darkgreen")
#lines(z28_aupr$rec, z28_aupr$prec, col="blue")
#legend("topright",lty=1,lwd=2,legend=c("AUPR=0.6955647","AUROC=0.4749878"))
#save.image('~/data/htseqcount.RData')
