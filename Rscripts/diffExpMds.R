##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Aug 2011 nanomed project (Sonia Human Th17 CARs)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

setwd("/Users/xichen/")

# source required functions
source("code/microarray_r_scripts/util.R")
library(samr)

# set paths
path.input <- "data/Microarray_1/output/"
path.output <- "data/Microarray_1/results/"

# set file names
file.data <- paste(sep="",path.input,"sonia_28_icos_bb_time_series_data_matrix_Dec_11_2012.RData")
file.tfs <- paste(sep="",path.input,"humanTFs.RData")

# load data (variable name is d)
load(file.data)
# load human tf names humanTFNames
load(file.tfs)

# get names of experiments we want to compare to ICOS 4hr to 0hr
icos.ctrl.ix <- grep("30h",colnames(d),value=T)## ICOS 0hr 
bb.ctrl.ix <- grep("20h",colnames(d),value=T)## BBz 0hr
z28.ctrl.ix <- grep("10h",colnames(d),value=T)## 28z 0hr
#h24.expt.ix <- grep("324h",colnames(d),value=T)## ICOS 24hr
#d4.expt.ix <- grep("34days",colnames(d),value=T)## ICOS 4days

## cd4.expt.names <- paste(sep="","cd4_t",c(2,6,12,24,48),"hr")
## cd8.expt.names <- paste(sep="","cd8_t",c(2,6,12,24,48),"hr")
icos.expt.names <- paste(sep="","3",c("4h","8h","24h","4days"))
bb.expt.names <- paste(sep="","2",c("4h","8h","24h","4days"))
z28.expt.names <- paste(sep="","1",c("4h","8h","24h","4days"))

all.expt.names <- c(icos.expt.names[1],bb.expt.names[1],z28.expt.names[1])
# create a list to hold diff expression results (and fold change results) of each comparison
#res <- list()
#fc <- list()

## data matrix to store p-value, t-test score, fold-change of differential expression between CARs  
w.p <- matrix(0,nr=dim(d)[1],nc=3)
rownames(w.p) <- rownames(d)
colnames(w.p) <- c(paste(sep="",icos.expt.names[1], "_vs_", bb.expt.names[1]),paste(sep="",bb.expt.names[1], "_vs_", z28.expt.names[1]),paste(sep="",icos.expt.names[1], "_vs_", z28.expt.names[1]))
  
w.t <- w.p
w.fc <- w.p

# create SAM data object (to compare icos expt_i to bbz control)
#for(i in 1:length(all.expt.names )){
  #time.point <- strsplit(icos.expt.names[i],"3")[[1]][2]
  #h0.expt.ix <- grep("30h",colnames(d),value=T)##ICOS 0hr 
  expt.ix <- grep(icos.expt.names[1],colnames(d),value=T)
  ctrl.ix <- paste(sep="",c("A","B","C"),"24h")#(bb.expt.names[1],colnames(d),value=T)
  d.sub <- d[,c(ctrl.ix,expt.ix)]
  ##two-class paired
  data <- list(x=d.sub,y=c(-(1:length(ctrl.ix)),(1:length(expt.ix))),
               geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  ##two-class unpaired
  #data <- list(x=d.sub,y=c(rep(1,length(icos.ctrl.ix)),rep(2,length(expt.ix))),
  #             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  sam.mat <- run.sam(data)
  w.p[rownames(sam.mat),paste(sep="",icos.expt.names[1], "_vs_", bb.expt.names[1])] <- 2*pt(-abs(sam.mat[,"t_test"]),df = 2)
  w.t[rownames(sam.mat),paste(sep="",icos.expt.names[1], "_vs_", bb.expt.names[1])] <- sam.mat[,"t_test"]
  w.fc[rownames(sam.mat),paste(sep="",icos.expt.names[1], "_vs_", bb.expt.names[1])] <- sam.mat[,"fold_change"]
  #fc[[ time.point ]][rownames(sam.mat),paste(sep="", icos.expt.names[i], "_vs_ctrl")] <- sam.mat[,"fold_change"]
#}

# create SAM data object (to compare icos expt_i to z28 control)
expt.ix <- grep(icos.expt.names[1],colnames(d),value=T)
ctrl.ix <- grep(z28.expt.names[1],colnames(d),value=T)
d.sub <- d[,c(ctrl.ix,expt.ix)]
##two-class paired
data <- list(x=d.sub,y=c(-(1:length(ctrl.ix)),(1:length(expt.ix))),
             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
##two-class unpaired
#data <- list(x=d.sub,y=c(rep(1,length(icos.ctrl.ix)),rep(2,length(expt.ix))),
#             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
sam.mat <- run.sam(data)
w.p[rownames(sam.mat),paste(sep="",icos.expt.names[1], "_vs_", z28.expt.names[1])] <- 2*pt(-abs(sam.mat[,"t_test"]),df = 2)
w.t[rownames(sam.mat),paste(sep="",icos.expt.names[1], "_vs_", z28.expt.names[1])] <- sam.mat[,"t_test"]
w.fc[rownames(sam.mat),paste(sep="",icos.expt.names[1], "_vs_", z28.expt.names[1])] <- sam.mat[,"fold_change"]

# create SAM data object (to compare bb expt_i to z28 control)
expt.ix <- paste(sep="",c("A","B","C"),"24h")
ctrl.ix <- grep(z28.expt.names[1],colnames(d),value=T)
d.sub <- d[,c(ctrl.ix,expt.ix)]
##two-class paired
data <- list(x=d.sub,y=c(-(1:length(ctrl.ix)),(1:length(expt.ix))),
             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
##two-class unpaired
#data <- list(x=d.sub,y=c(rep(1,length(icos.ctrl.ix)),rep(2,length(expt.ix))),
#             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
sam.mat <- run.sam(data)
w.p[rownames(sam.mat),paste(sep="",bb.expt.names[1], "_vs_", z28.expt.names[1])] <- 2*pt(-abs(sam.mat[,"t_test"]),df = 2)
w.t[rownames(sam.mat),paste(sep="",bb.expt.names[1], "_vs_", z28.expt.names[1])] <- sam.mat[,"t_test"]
w.fc[rownames(sam.mat),paste(sep="",bb.expt.names[1], "_vs_", z28.expt.names[1])] <- sam.mat[,"fold_change"]

j <- 0
for(i in 1:length(rownames(w.p))) {
	if(w.p[i,1] < 0.05 | w.p[i,2] < 0.05 | w.p[i,3] < 0.05){
		j = j+1
	}
}

n.icos.bbz <- 0
for(i in 1:length(rownames(w.p))) {
	if(w.p[i,1] < 0.05 & abs(log2(w.fc[i,1])) > 1){
		n.icos.bbz <- n.icos.bbz +1
	}
}

n.bbz.z28 <- 0
for(i in 1:length(rownames(w.p))) {
	if(w.p[i,2] < 0.05  &  abs(log2(w.fc[i,2])) > 1){
		n.bbz.z28 <- n.bbz.z28 +1
	}
}

n.icos.z28 <- 0
for(i in 1:length(rownames(w.p))) {
	if(w.p[i,3] < 0.05 & abs(log2(w.fc[i,3])) > 1){
		n.icos.z28 <- n.icos.z28 +1
	}
}

res <- matrix(0,nr=j,nc=3)
colnames(res) <- colnames(w.p)
k <- 1
for(i in 1:length(rownames(w.p))) {
	if(w.p[i,1] < 0.05 | w.p[i,2] < 0.05 | w.p[i,3] < 0.05) {
		res[k,] <- w.p[i,]
		if (k == 1){
			n <- rownames(w.p)[i]
		}else{
		n <- c(n,rownames(w.p)[i])
		}
		k = k+1
	}
}

rownames(res) <- n

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
x <- read.delim("~/data/Microarray_1/output/sonia_28_icos_bb_time_series_data_matrix_Dec_11_2012.xls", header=T)
x <- as.matrix(x)
# grep for all cars all donors only 4hr
ix <- grep("^..4h",colnames(x),perl=T,value=T)

x.sub <- x[n,]
x.sub.sub <- x.sub[,ix]
x.sub.sub.t <- t(x.sub.sub)
d <- dist(x.sub.sub.t) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
nms.vec <- gsub("4h","",rownames(fit$points))
nms.vec <- gsub("2","-BBz",nms.vec)
nms.vec <- gsub("1","-28z",nms.vec)
nms.vec <- gsub("3","-ICOSz",nms.vec)

library("vegan")
x <- fit$points[,1]
y <- fit$points[,2]
plot(fit$points,col=c("gray","black","red","gray","black","red","gray","black","red"),pch=20,xlab = "Dimension 1",ylab = "Dimension 2",xlim = c(-20,20),ylim = c(-20,20),lwd = 5)
legend("topright",pch = 20,lwd=2,legend=c("ICOSz","BBz","28z"),col=c("red","black","gray"),bty = "n")
text(x, y, labels = c("A","A","A","B","B","B","C","C","C"), cex=1,pos=3)
# ordiellipse(fit$points[c(3,6,9),],c(1,1,1),conf=0.9, kind = "sd",lwd=2, draw = "polygon", border = "red")
ordiellipse(fit$points[c(2,5,8),],c(1,1,1),conf=0.9, kind = "sd",lwd=2, draw = "polygon", border = "black")
ordiellipse(fit$points[c(1,4,7),],c(1,1,1),conf=0.9, kind = "sd",lwd=2, draw = "polygon", border = "gray")

xy.icos <- fit$points[c(3,6,9),]
exy.icos <- ellipsoidhull(as.matrix(xy.icos))
lines(predict(exy.icos),col = "red",lwd = 4)

exy.icos <- ellipsoidhull(as.matrix(xy.icos))
exy.icos$cov[1,1] <- exy.icos$cov[1,1]*1.5
exy.icos$cov[2,2] <- exy.icos$cov[2,2]*3
lines(predict(exy.icos),col = "red",lwd = 2)

xy.z28 <- fit$points[c(1,4,7),]
exy.z28 <- ellipsoidhull(as.matrix(xy.z28))
lines(predict(exy.z28),col = "gray",lwd = 4)

xy.bbz <- fit$points[c(2,5,8),]
exy.bbz <- ellipsoidhull(as.matrix(xy.bbz))
lines(predict(exy.bbz),col = "black",lwd = 4)

d.icos.bbz <- sqrt((exy.icos$loc[1]-exy.bbz$loc[1])^2 + (exy.icos$loc[2]-exy.bbz$loc[2])^2)
d.icos.z28 <- sqrt((exy.icos$loc[1]-exy.z28$loc[1])^2 + (exy.icos$loc[2]-exy.z28$loc[2])^2)
d.bbz.z28 <- sqrt((exy.bbz$loc[1]-exy.z28$loc[1])^2 + (exy.bbz$loc[2]-exy.z28$loc[2])^2)
#######################################^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



# calc specificity of diff exression between different CARs
for(i in 1:length(icos.expt.names )){
  time.point <- strsplit(icos.expt.names[i],"3")[[1]][2]
  ## calc difference btwn icos and bb/z28 lists of differentially expressed genes (for 48hr vs. ctrl)
  res[[ time.point ]][,"diff_icos_vs_bb"] <- res[[ time.point ]][,1]-res[[ time.point ]][,2]
  fc[[ time.point ]][,"diff_icos_vs_bb"] <- fc[[ time.point ]][,1]-fc[[ time.point ]][,2]
  res[[ time.point ]][,"diff_icos_vs_z28"] <- res[[ time.point ]][,1]-res[[ time.point ]][,3]
  fc[[ time.point ]][,"diff_icos_vs_z28"] <- fc[[ time.point ]][,1]-fc[[ time.point ]][,3]
  ## for each tf put a 1 for non tf put a 0 in column is.tf
  ix <- which(rownames(w) %in% humanTFNames)
  res[[ time.point ]][ix,"is_tf"] <- 1
  fc[[ time.point ]][ix,"is_tf"] <- 1
}
x <- res[["4days"]]
y <- fc[["4days"]]
write.table(x,sep="\t",,file=paste(sep="",path.input,"sam_diff_exp_4days_vs_ctrl_paired.xls"))
write.table(y,sep="\t",,file=paste(sep="",path.input,"fold_change_4days_vs_ctrl_paired.xls"))

# get a matrix with specificiy scores at each time point (positive means upregulated in cd4, negative means upregulated in cd8)
w <- matrix(0,nr=dim(res[[ 1 ]])[1],nc=length(res))
rownames(w) <- rownames(res[[ 1 ]])
colnames(w) <- names(res)
for(i in 1:length(res)){
  time.point <- names(res)[i]
  #w[,time.point] <- res[[time.point]][,"diff_icos_vs_bb"]
  w[,time.point] <- res[[time.point]][,"diff_icos_vs_z28"]
}
w <- w[order(abs(w[,"4h"]),decreasing=T),]

w.fc <- matrix(0,nr=dim(res[[ 1 ]])[1],nc=length(res))
rownames(w.fc) <- rownames(res[[ 1 ]])
colnames(w.fc) <- names(res)
for(i in 1:length(fc)){
  time.point <- names(res)[i]
  #w.fc[,time.point] <- fc[[time.point]][,"diff_icos_vs_bb"]
  w.fc[,time.point] <- fc[[time.point]][,"34h_vs_ctrl"]
  #w.fc[,time.point] <- fc[[time.point]][,"24h_vs_ctrl"]
  #w.fc[,time.point] <- fc[[time.point]][,"14h_vs_ctrl"]
  #w.fc[,time.point] <- fc[[time.point]][,"diff_icos_vs_z28"]
}
w.fc <- w.fc[order(abs(w.fc[,"4h"]),decreasing=T),]
#write.table(w.fc,sep="\t",,file=paste(sep="",path.input,"fold_change_icos_vs_bb_unpaired.xls"))
write.table(w.fc,sep="\t",,file=paste(sep="",path.input,"fold_change_icos_4h_vs_ctrl_paired.xls"))
#write.table(w.fc,sep="\t",,file=paste(sep="",path.input,"fold_change_bb_4h_vs_ctrl_paired.xls"))
#write.table(w.fc,sep="\t",,file=paste(sep="",path.input,"fold_change_z28_4h_vs_ctrl_paired.xls"))


cut.sam <- 4
cut.fc <- 2
gns.sam <- rownames(w)[which(w[,"4h"]>cut.sam)]
gns.fc <- rownames(w)[which(w.fc[,"4h"]>cut.fc)]
gns.volcano <- intersect(gns.sam,gns.fc)

# take genes with big difference in expression at time 48hr
expt.ix.icos<- grep(icos.expt.names[1],colnames(d),value=T)
expt.ix.bb<- grep(paste(sep="","^.",bb.expt.names[i]),colnames(d),value=T)##grep(bb.expt.names[1],colnames(d),value=T)
x=d[,expt.ix.icos]-d[,expt.ix.bb]
x.mean=apply(x,1,mean)
x.sum=apply(x,1,sum)
ix=sort(x.mean,decreasing=T,index.return=T)$ix
gns.abs.diff.4h <- names(x.mean)[ix]


gns.intersting <- c("IL23R","IL21","PRG4","AIM2","IL1R1","FNBP1L","IL1A","IL17F","NEK3","ETV6","FRMD4B","CXCR5","NCKAP1","PPP4R4","SPATS2L",
                    "CD40LG","CD4","RTKN2","FAM40B","GPR87","HSD11B1","CTLA4","CD109","PTGR1",
                    "IL1R2","CPM","LY75","CCL20","DYNC2LI1","IL2","MYOF","DEPDC1","PLEKHH2","SNORA73A","CYFIP1")
## gns <- gns.abs.diff.48hr[1:100]
## gns <- c("ETV6","SKIL","ZNF670","EGR1","AHR","IL23R","IL21","CCR2","TSPYL2","CCR7")
gns <- gns.intersting
pdf(file=paste(sep="",path.output,"interesting_genes.pdf"))
for(i in 1:length(gns)){
  plot.gn(d,w,gns[i],cex=1.5)
}
dev.off()


gns.intersting <- c("IL23R","IL21","PRG4","AIM2","IL1R1","FNBP1L","IL1A","IL17F","NEK3","ETV6","FRMD4B","CXCR5","NCKAP1","PPP4R4","SPATS2L",
                    "CD40LG","CD4","RTKN2","FAM40B","GPR87","HSD11B1","CTLA4","CD109","PTGR1",
                    "IL1R2","CPM","LY75","CCL20","DYNC2LI1","IL2","MYOF","DEPDC1","PLEKHH2","SNORA73A","CYFIP1")