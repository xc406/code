##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Aug 2011 nanomed project (Mike Milones human cd4/cd8 data)
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

# remove three bad experiments pt3_cd8_t0hr,pt3_cd8_t6hr, pt3_cd8_t24hr
## ix.rm <- which(colnames(d) %in% c("pt3_cd8_t0hr", "pt3_cd8_t6hr", "pt3_cd8_t24hr"))
## d <- d[,-ix.rm]

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

all.expt.names <- c(icos.expt.names,bb.expt.names,z28.expt.names)
# create a list to hold diff expression results (and fold change results) of each comparison
res <- list()
fc <- list()
for(i in 1:length(icos.expt.names)){
  time.point <- strsplit(icos.expt.names[i],"3")[[1]][2]
  w <- matrix(0,nr=dim(d)[1],nc=6)
  rownames(w) <- rownames(d)
  colnames(w) <- c(paste(sep="",icos.expt.names[i], "_vs_ctrl"),paste(sep="",bb.expt.names[i], "_vs_ctrl"),paste(sep="",z28.expt.names[i], "_vs_ctrl"),"diff_icos_vs_bb","diff_icos_vs_z28","is_tf")
  res[[ time.point ]] <- w
  fc[[ time.point ]] <- w
}

# create SAM data object (to compare icos expt_i to its control)
for(i in 1:length(icos.expt.names )){
  time.point <- strsplit(icos.expt.names[i],"3")[[1]][2]
  #h0.expt.ix <- grep("30h",colnames(d),value=T)##ICOS 0hr 
  expt.ix <- grep(icos.expt.names[i],colnames(d),value=T)
  
  d.sub <- d[,c(icos.ctrl.ix,expt.ix)]
  ##two-class paired
  data <- list(x=d.sub,y=c(-(1:length(icos.ctrl.ix)),(1:length(expt.ix))),
               geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  ##two-class unpaired
  #data <- list(x=d.sub,y=c(rep(1,length(icos.ctrl.ix)),rep(2,length(expt.ix))),
  #             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  sam.mat <- run.sam(data)
  res[[ time.point ]][rownames(sam.mat),paste(sep="",icos.expt.names[i], "_vs_ctrl" )] <- sam.mat[,"t_test"]
  fc[[ time.point ]][rownames(sam.mat),paste(sep="", icos.expt.names[i], "_vs_ctrl")] <- sam.mat[,"fold_change"]
}

# create SAM data object (to compare bbz expt_i to its control)
for(i in 1:length(bb.expt.names )){
  time.point <- sub('^2', '', bb.expt.names[i])##strsplit(bb.expt.names[i],"^2")[[1]][2] ##substitute first char with nothing!
  expt.ix <- grep(paste(sep="","^.",bb.expt.names[i]),colnames(d),value=T)
  d.sub <- d[,c(bb.ctrl.ix,expt.ix)]
  #two-class paired 
  data <- list(x=d.sub,y=c(-(1:length(bb.ctrl.ix)),(1:length(expt.ix))),
               geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  #two-class unpaired
  #data <- list(x=d.sub,y=c(rep(1,length(bb.ctrl.ix)),rep(2,length(expt.ix))),
  #             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  sam.mat <- run.sam(data)
  res[[ time.point ]][rownames(sam.mat),paste(sep="",bb.expt.names[i], "_vs_ctrl")] <- sam.mat[,"t_test"]
  fc[[ time.point ]][rownames(sam.mat),paste(sep="",bb.expt.names[i], "_vs_ctrl")] <- sam.mat[,"fold_change"]
}

# create SAM data object (to compare 28z expt_i to its control)
for(i in 1:length(z28.expt.names )){
  time.point <- strsplit(z28.expt.names[i],"1")[[1]][2]
  #h0.expt.ix <- grep("30h",colnames(d),value=T)##ICOS 0hr 
  expt.ix <- grep(z28.expt.names[i],colnames(d),value=T)
  
  d.sub <- d[,c(z28.ctrl.ix,expt.ix)]
  ##two-class paired
  data <- list(x=d.sub,y=c(-(1:length(z28.ctrl.ix)),1:length(expt.ix)),
               geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  ##two-class unpaired
  #data <- list(x=d.sub,y=c(rep(1,length(z28.ctrl.ix)),rep(2,length(expt.ix))),
  #             geneid=as.character(1:dim(d)[1]),genenames=rownames(d), logged2=TRUE)
  sam.mat <- run.sam(data)
  res[[ time.point ]][rownames(sam.mat),paste(sep="",z28.expt.names[i], "_vs_ctrl" )] <- sam.mat[,"t_test"]
  fc[[ time.point ]][rownames(sam.mat),paste(sep="", z28.expt.names[i], "_vs_ctrl")] <- sam.mat[,"fold_change"]
}


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
x <- res[["4h"]]
y <- fc[["4h"]]
write.table(x,sep="\t",,file=paste(sep="",path.input,"sam_diff_exp_4h_vs_ctrl_paired.xls"))
write.table(y,sep="\t",,file=paste(sep="",path.input,"fold_change_4h_vs_ctrl_paired.xls"))







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
