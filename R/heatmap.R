library(pheatmap)
path.input <- "/Users/xclair/Documents/google-python-exercises/"

# color scheme:
or.col <- c(255,102,0)
bl.col <- c(0,51,255)
blue.col <- rgb(bl.col[1],bl.col[2],bl.col[3],max=255)
orange.col <- rgb(or.col[1],or.col[2],or.col[3],max=255)

f.nm <- paste(sep="",path.input,"Signature_list_for_heatmaphg.txt")
x <- read.table(f.nm,sep="\t",header=F,as.is=T)
colnames(x) <- c("","Th1","Th2","Th17","Tfh","iTreg")

# put fold change values here
x.fc <- read.delim(paste(sep="",path.input,"fold_change_4h_vs_ctrl_paired.xls"), header=T)
x.fc <- as.matrix(x.fc)
# now extract icos = "X34h_vs_ctrl"
fc.icos <- log2(x.fc[,"X34h_vs_ctrl"])##icos 4h vs 0h 
# now extract bbz = "X24h_vs_ctrl"
fc.bbz <- log2(x.fc[,"X24h_vs_ctrl"])##bbz 4h vs 0h 
# now extract z28 = "X14h_vs_ctrl"
fc.z28 <- log2(x.fc[,"X14h_vs_ctrl"])##z28 4h vs 0h 


# determine boundaries (not so important)
#eff.max.val <- min( abs( c( max(x.fc), min(x.fc) ) ) )
eff.max.val <- 2
ix.icos <- which(fc.icos < -eff.max.val)
if(length(ix.icos)){fc.icos[ix.icos] <- -eff.max.val}
ix.icos <- which(fc.icos > eff.max.val)
if(length(ix.icos)){fc.icos[ix.icos] <- eff.max.val}

ix.bbz <- which(fc.bbz < -eff.max.val)
if(length(ix.bbz)){fc.bbz[ix.bbz] <- -eff.max.val}
ix.bbz <- which(fc.bbz > eff.max.val)
if(length(ix.bbz)){fc.bbz[ix.bbz] <- eff.max.val}

ix.z28 <- which(fc.z28 < -eff.max.val)
if(length(ix.z28)){fc.z28[ix.z28] <- -eff.max.val}
ix.z28 <- which(fc.z28 > eff.max.val)
if(length(ix.z28)){fc.z28[ix.z28] <- eff.max.val}

fc <- cbind(fc.icos,fc.bbz,fc.z28)

#fl.nm <- "your file name here"
#pdf(fl.nm)
#m <- m.cut # here goes the fc data
for(j in 2:ncol(x)){
	nm <- colnames(x)[j]
	gns <- x[,1][which(x[,j]==1)]
	
	# cluster the rows
	d <- dist(as.matrix(fc[gns,]))
	hc <- hclust(d)
	gn.order <- hc$labels[hc$order]
	
# plot heatmap
	if(length(gns)>1){
		pheatmap(fc[gn.order,],scale="none",color=colorRampPalette(c(orange.col,"white",
		blue.col))(255),breaks=seq(-1*eff.max.val,eff.max.val,by=(2*eff.max.val)/255),cellwidth = 36, cellheight = 12,treeheight_col=0,
		cex=1,cluster_cols=F,cluster_rows=F,main=colnames(x)[j])
					
	}
	
#	pheatmap(x.fc[gn.order],scale="none",cellheight = 4,main="FC",breaks=seq(-eff.max.val,eff.max.val,by=(2*eff.max.val)/tot ),
#		colorRampPalette(c(orange.col,"white",blue.col))(tot),	
#	cex=0.3,main=paste(sep="","kc=",cut.abs," pval=",deseq.pval.cut),cluster_rows=F,cluster_cols=F)
	
}
#dev.off()