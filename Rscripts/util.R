##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Aug 2011 nanomed project (Mike Milones human cd4/cd8 data)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

plot.gn <- function(d,s.mat,gn.nm,cex=2){
  par(mfrow=c(2,1))

  ## plot time series for cd4cex
  x1 <- d[gn.nm,grep("pt1_cd4",colnames(d),perl=T)]
  x2 <- d[gn.nm,grep("pt2_cd4",colnames(d),perl=T)]
  x3 <- d[gn.nm,grep("pt3_cd4",colnames(d),perl=T)]
  y.max <- max(c(x1,x2,x3))
  y.min <- min(c(x1,x2,x3))
  plot(x1,pch=20,col="darkgreen",ylim=c(y.min,y.max),xaxt="n",type="l",xlab="",ylab="cd4",main=gn.nm,cex.main=cex,cex.lab=cex,cex.axis=cex)
  points(x2,pch=20,col="darkblue",type="l")
  points(x3,pch=20,col="darkred",type="l")
  axis(1, at=1:length(x1),labels=sapply(strsplit(names(x1),"_"),function(i) i[3]), las=2,cex.axis=cex)
  legend("topleft",inset=.01,title="patient",c("pt1","pt2","pt3"),fill=c("darkgreen","darkblue","darkred"))

  ## plot time series for cd4
  x1 <- d[gn.nm,grep("pt1_cd8",colnames(d),perl=T)]
  x2 <- d[gn.nm,grep("pt2_cd8",colnames(d),perl=T)]
  x3 <- d[gn.nm,grep("pt3_cd8",colnames(d),perl=T)]
  y.max <- max(c(x1,x2,x3))
  y.min <- min(c(x1,x2,x3))
  plot(x1,pch=20,col="darkgreen",ylim=c(y.min,y.max),xaxt="n",type="l",xlab="",ylab="cd8",cex.axis=cex,cex.lab=cex)
  points(x2,pch=20,col="darkblue",type="l")
  points(x3,pch=20,col="darkred",type="l")
  axis(1, at=1:length(x1),labels=sapply(strsplit(names(x1),"_"),function(i) i[3]), las=2,cex.axis=cex)
  legend("topleft",inset=.01,title="patient",c("pt1","pt2","pt3"),fill=c("darkgreen","darkblue","darkred"))

  ## barplot(s.mat[gn.nm,],ylab="specificity [z(cd4)-z(cd8)]",cex.axis=cex,cex.lab=cex,cex=cex)
}

## function runs SAM on object of type data
# input:
#       - data: list of:
#         - dataset x
#         - vector of experiment assigment (e.g. 1,1,1,2,2,2)
#         - vector of gene ids = 1,2,...,N
#         - vector of gene names = gn.name1,gn.name2,...,genenameN
#         - logical stating if data was log2 transformed (TRUE)
# output:
#      - sam.mat: a matrix with two columns (zscore,and q value) and a row for each gene

run.sam <- function(data){
  samr.obj <- samr(data,resp.type="Two class paired",nperms=100)
  delta=-1
  samr.plot(samr.obj,delta)
  delta.table <- samr.compute.delta.table(samr.obj)
  siggenes.table <-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)
  dval.ranks.mat <- rbind(as.matrix(as.numeric(siggenes.table$genes.up[,"Score(d)"])),as.matrix(as.numeric(siggenes.table$genes.lo[,"Score(d)"])))
  qval.ranks.mat <- rbind(as.matrix(as.numeric(siggenes.table$genes.up[,"q-value(%)"])),as.matrix(as.numeric(siggenes.table$genes.lo[,"q-value(%)"])))
  fc.ranks.mat <- rbind(as.matrix(as.numeric(siggenes.table$genes.up[,"Fold Change"])),as.matrix(as.numeric(siggenes.table$genes.lo[,"Fold Change"])))
  sam.mat <- cbind(dval.ranks.mat,fc.ranks.mat,qval.ranks.mat)
  colnames(sam.mat) <- c("t_test","fold_change","q_value")
  rownames(sam.mat) <- c(siggenes.table$genes.up[,"Gene ID"],siggenes.table$genes.lo[,"Gene ID"])
  return(sam.mat)
}
