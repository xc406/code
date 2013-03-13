rm(list=ls())
## gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"

#load("/Users/xichen/data/results_combine_data_exprsn_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1.Rdata")
#pred.file <- knockout_rna_zscores.all <- read.delim(gs.file)

#k_r_network_stat3 <- res[["activation"]][["STAT3"]][["th17"]][,"KCRI"]+res[["activation"]][["BATF"]][["th17"]][,"KCRI"]+res[["activation"]][["IRF4"]][["th17"]][,"KCRI"]+res[["activation"]][["MAF"]][["th17"]][,"KCRI"]+res[["activation"]][["RORC"]][["th17"]][,"KCRI"]
#write.table(k_r_network_stat3,sep="\t",,file=paste(sep="","~/data/Microarray_1/output/kcri_core.xls"))
#k_r_network_stat3 <- k_r_network_stat3[which(k_r_network_stat3>0)]

gs <- read.delim("~/data/Microarray_1/output/gshg.txt", header=F)##74 genes as literature gs in mm9

#k_r_network_stat3 <- read.delim("~/data/Microarray_1/output/kcri_corehg.txt", header=F)
#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/Microarray_1/output/fold_change_4h_vs_ctrl_unpairedmm.txt", header=F)

# make sure with str, you get what you expect
motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/Microarray_1/output/sam_diff_exp_4h_vs_ctrl_paired.xls", header=T)
motif.and.dnase.reads.stat3.rep1 <- as.matrix(motif.and.dnase.reads.stat3.rep1)
# for consistancy lets have gene names in upper case form
rownames(motif.and.dnase.reads.stat3.rep1) <- toupper(rownames(motif.and.dnase.reads.stat3.rep1))

#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/Microarray_1/output/fold_change_4h_vs_ctrl_pairedmm.txt", header=F)

#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/htseq_output/final/Irf4_htseq_output_avg", header=F)
#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/phylopStat3Window_htseq_output", header=F)
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[-(nrow(motif.and.dnase.reads.stat3.rep1) - 0:4),] #remove meta
#motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>300),2] <- 300
#motif.and.dnase.reads.stat3.rep1 <- motif.and.dnase.reads.stat3.rep1[which(motif.and.dnase.reads.stat3.rep1[,2]>2.0), ]

# now extract icos = "X34h_vs_ctrl"
sam_diff_exp_icos <- motif.and.dnase.reads.stat3.rep1[,"X34h_vs_ctrl"]##icos 4h vs 0h 
# now extract bbz = "X24h_vs_ctrl"
sam_diff_exp_bbz <- motif.and.dnase.reads.stat3.rep1[,"X24h_vs_ctrl"]##icos 4h vs 0h 
# now extract z28 = "X14h_vs_ctrl"
sam_diff_exp_z28 <- motif.and.dnase.reads.stat3.rep1[,"X14h_vs_ctrl"]##icos 4h vs 0h 


#sam_diff_exp_icos <- motif.and.dnase.reads.stat3.rep1[,2]##icos 4h vs 0h 
#sam_diff_exp_bbz <- motif.and.dnase.reads.stat3.rep1[,3]##bbz 4h vs 0h
#sam_diff_exp_z28 <- motif.and.dnase.reads.stat3.rep1[,4]##28z 4h vs 0h 
#names(sam_diff_exp_icos) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])
#names(sam_diff_exp_bbz) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])
#names(sam_diff_exp_z28) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])

#k_r <- k_r_network_stat3[,2]##28z 4h vs 0h 
#names(k_r) <- toupper(k_r_network_stat3[,1])
#m_d_network <- sort(k_r, decreasing = T)
#m_d <- c(rep(NA,2000))

#m_d <- m_d_network[1:2000]
#k_r_network_stat3 <- order(k_r_network_stat3, decreasing = T)

#m_d_network_stat3 <- order(m_d_network_stat3, decreasing = T)
#m_d_network_stat3 <- motif.and.dnase.reads.stat3.rep1[,2]

th17_network_stat3<-c(rep(1,length(k_r_network_stat3[,1])))
names(th17_network_stat3) <- toupper(k_r_network_stat3[,1])

#unique.names <- unique(c( names(m_d_network_stat3), names(k_r_network_stat3)))
#unique.names <- unique(c( names(m_d), names(k_r_network_stat3)))
#unique.names <- unique(c( names(m_d), names(sam_diff_exp_icos)))
#unique.names <- unique(c( names(sam_diff_exp_icos), names(th17_network_stat3)))
gs.genes <- intersect(c( names(sam_diff_exp_icos), names(th17_network_stat3)))


#overlap <- (length(m_d_network_stat3) + length(k_r_network_stat3)) - length(unique.names)
data <- matrix(NA, nr=length(unique.names), nc=4)
rownames(data) <- unique.names
colnames(data) <- c("gs_74genes", "sam_diff_exp_icos","sam_diff_exp_bbz","sam_diff_exp_z28")
for (i in 1:nrow(data)) {
  
  data[i,2] <- sam_diff_exp_icos[unique.names[i]]##icos sam_diff_exp 4h vs 0h
  data[i,3] <- sam_diff_exp_bbz[unique.names[i]]##bbz 
  data[i,4] <- sam_diff_exp_z28[unique.names[i]]##28z 
  #data[i,1] <- m_d[unique.names[i]]##k_r_network_stat3[unique.names[i]]##kcri network score
  data[i,1] <- th17_network_stat3[unique.names[i]]##74 th17 genes from literature
}
#data[which(is.na(data))] <- 0

gs <- data[,1]
gs[which(is.na(gs))] <- 0 
gs <- gs>0

pred.icos <- data[,2]
pred.bbz <- data[,3]
pred.z28 <- data[,4]
pred.icos <- pred.icos[which(!is.na(pred.icos))]
pred.bbz <- pred.bbz[which(!is.na(pred.bbz))]
pred.z28 <- pred.z28[which(!is.na(pred.z28))]

############ Aviv touchups ##########
pred.icos.sort <- sort(pred.icos,decreasing=T)
pred.bbz.sort <- sort(pred.bbz,decreasing=T)
pred.z28.sort <- sort(pred.z28,decreasing=T)
#get vecs for gsea
gene.list.icos <- names(pred.icos.sort) # names
gene.list.bbz <- names(pred.bbz.sort)
gene.list.z28 <- names(pred.z28.sort)
correl.vec.icos <- pred.icos.sort # values
correl.vec.bbz <- pred.bbz.sort
correl.vec.z28 <- pred.z28.sort

gene.set <- names(gs)[which(gs == 1)]
source('~/code/GSEA.R')
GSEA.res.icos <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.icos, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.icos)
GSEA.res.bbz <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.bbz, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.bbz)
GSEA.res.z28 <- GSEA.EnrichmentScore(run.perm=FALSE,gene.list=gene.list.z28, gene.set=gene.set,weighted.score.type = 1, correl.vector = correl.vec.z28)

source('~/code/AUPR/aupr.R')
AUPR.res.icos <- calcAupr(pred.icos,gs)
AUPR.res.bbz <- calcAupr(pred.bbz,gs)
AUPR.res.z28 <- calcAupr(pred.z28,gs)
############ Aviv END ##########

##AUPR plot
plot.new()
plot(AUPR.res.icos$rec, AUPR.res.icos$prec, type = 'lines', main = 'CARs_4h_vs_0h_diff_exp vs 74_Th17_genes', xlab = 'recall', ylab = 'precision',xlim = c(0,0.3), col="red")
#plot(AUPR.res.icos$rec, AUPR.res.icos$prec, type = 'lines', main = 'CARs_4h_vs_0h_diff_exp vs KCRI network', xlab = 'recall', ylab = 'precision',xlim = c(0,0.3), col="red")
lines(AUPR.res.bbz$rec, AUPR.res.bbz$prec, col="darkgreen")
lines(AUPR.res.z28$rec, AUPR.res.z28$prec, col="blue")
#legend("topright",lty=1,lwd=2,legend=c("AUPR.icos=0.03649203","AUPR.bbz=0.02804974","AUPR.28z=0.02672071"),col=c("red","darkgreen","blue"))
legend("topright",lty=1,lwd=2,legend=c(paste(sep="","AUPR.icos=", AUPR.res.icos$AUPR),paste(sep="","AUPR.bbz=", AUPR.res.bbz$AUPR),paste(sep="","AUPR.z28=", AUPR.res.z28$AUPR)),col=c("red","darkgreen","blue"))
##AUROC plot
plot.new()
#plot(AUPR.res.icos$fpr, AUPR.res.icos$rec, type = 'lines', main = 'CARs_4h_vs_0h_diff_exp vs 74_Th17_genes', xlab = 'false positive rate', ylab = 'true positive rate', col="red")
plot(AUPR.res.icos$fpr, AUPR.res.icos$rec, type = 'lines', main = 'CARs_4h_vs_0h_diff_exp vs KCRI network', xlab = 'false positive rate', ylab = 'true positive rate', col="red")
lines(AUPR.res.bbz$fpr, AUPR.res.bbz$rec, col="darkgreen")
lines(AUPR.res.z28$fpr, AUPR.res.z28$rec, col="blue")
#legend("topleft",lty=1,lwd=2,legend=c("AUROC.icos=0.7015295","AUROC.bbz=0.6703645","AUROC.28z=0.6561663"),col=c("red","darkgreen","blue"))
#legend("topleft",lty=1,lwd=2,legend=c("AUROC.icos","AUROC.bbz","AUROC.28z"),col=c("red","darkgreen","blue"))
legend("topright",lty=1,lwd=2,legend=c(paste(sep="","AUROC.icos=", AUPR.res.icos$AUROC),paste(sep="","AUROC.bbz=", AUPR.res.bbz$AUROC),paste(sep="","AUROC.z28=", AUPR.res.z28$AUROC)),col=c("red","darkgreen","blue"))
# gs <- data[,2]
# colnames(data)
# pred <- data[,1]
# gs>0
# pred <- pred>0
# gs
#calcAupr(pred,gs)
#icosp_aupr <- calcAupr(pred,gs)
#plot(icosp_aupr$rec, icosp_aupr$prec, type = 'lines', main = 'ICOSz_4h_vs_0h_vs_74genes', xlab = 'recall', ylab = 'precision')
#lines(bb_aupr$rec, bb_aupr$prec, col="darkgreen")
#lines(z28_aupr$rec, z28_aupr$prec, col="blue")
#legend("topright",lty=1,lwd=2,legend=c("icos","bbz","28z"),col=c("darkred","darkgreen","blue"))
#save.image('~/data/htseqcount.RData')

res <- list(GSEA.res.icos,GSEA.res.bbz,GSEA.res.z28)
names(res) <- c("icos","bbz","z28")
#f.nm <- paste(sep="",path.output,"gsea_",gs.lists[iter],"_",comb.case,".pdf")
#pdf(f.nm)
cls <- c("firebrick3",rgb(255, 106, 106,maxColorValue=255, 150),rgb(154, 205, 50,maxColorValue=255, 150),"darkgray")
nf <- layout(rbind(1,2,3),heights=c(3,1,2))
op <- par(mar=c(0,7,0.5,0.5),mgp=c(4, 1, 0))

for(i in 1:length(res)){
	n <- length(res[[i]]$RES)
	y <- res[[i]]$RES
	x <- 1:n
	if(i==1){
		plot(x=x,y=y,type="l",lwd=3,col=cls[i],cex.lab=2,ylim=c(-0.2,1),ylab="Enrichment Score (ES)",xaxt="n",xlab="",cex.axis=2,las=1)
	} else {
		lines(x=x,y=y,lwd=3,col=cls[i],cex.axis=2,las=1)
	}
}
legend(y=1,x=(n-4000),legend=c("ICOS","BBz","z28"),lty=1,lwd=3,col=cls[1:length(res)],cex=2)

coords <- list()
coords[[3]] <- c(-1,-0.5)
coords[[2]] <- c(-0.25,0.25)
coords[[1]] <- c(0.5,1)
for(k in 1:length(res)){
	if(k ==1){
		plot(y=c(-0.05,0.05),x=c(1,1),col="gray89",type="l",lty=1,xlim=c(0,n),lwd=0.1, xlab="",ylab="",axes=F,ylim=c(-1,1))
	}
	x <- c(0,0,n,n)
	y <- c(coords[[k]],rev(coords[[k]]))
	polygon(x,y, col="gray89")
	s <- res[[k]]$indicator
	ix <- which(s == 1)
	for(i in ix){
		lines(y=coords[[k]],x=c(i,i),col="darkred",lwd=2)
	}
}
for(k in 1:length(res)){
	i = res[[k]]$arg.ES
	y <- c(coords[[k]])
	lines(y=y,x=c(i,i),col="darkblue",lwd=5)
}

mtext("ICOS",side=2,cex=1,las=1,line = -1,at=0.75)
mtext("BBz",side=2,cex=1,las=1,line = -1,at=0)
mtext("z28",side=2,cex=1,las=1,line = -1,at=-0.75)

par(mar=c(4.5,7,0.5,0.5),mgp=c(3, 1, 0))
de <- sort(sam_diff_exp_icos,decreasing=T)
plot(y=de,x = c(1:length(res[[1]]$RES)),col=cls[1],pch=20,lty=1,xlim=c(0,n),cex.axis=1,cex.lab=1.2,las=1,lwd=0.5,ylab="FC",xlab="Ranked gene list",ylim=c(min(de),max(de)))
legend(y=max(de)-1,x=(n-8000),legend=c("ICOS","bbz","z28"),lty=1,lwd=3,col=cls[c(1,2,3)],cex=1)
lines(x=1:n,y=sort(sam_diff_exp_bbz,decreasing=T),pch=20,col=cls[2],lwd=1)
lines(x=1:n,y=sort(sam_diff_exp_z28,decreasing=T),pch=20,col=cls[3],lwd=1)
