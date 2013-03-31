
path.input <- "/Users/xichen/data/Microarray_1/output/"

# color scheme:
or.col <- c(255,102,0)
bl.col <- c(0,51,255)
blue.col <- rgb(bl.col[1],bl.col[2],bl.col[3],max=255)
orange.col <- rgb(or.col[1],or.col[2],or.col[3],max=255)

f.nm <- paste(sep="",path.input,"Signature_list_sonia_heatmaphg.txt")
x <- read.table(f.nm,sep="\t",header=F,as.is=T)
colnames(x) <- c("","Th1","Th2","Th17","iTreg")

# put fold change values here
x.fc <- read.delim(paste(sep="",path.input,"fold_change_4h_vs_ctrl_paired.xls"), header=T)
x.fc <- as.matrix(x.fc)
# now extract icos = "X34h_vs_ctrl"
fc.icos <- log2(as.numeric(x.fc[,"X34h_vs_ctrl"]))##icos 4h vs 0h 
# now extract bbz = "X24h_vs_ctrl"
fc.bbz <- log2(as.numeric(x.fc[,"X24h_vs_ctrl"]))##bbz 4h vs 0h 
# now extract z28 = "X14h_vs_ctrl"
fc.z28 <- log2(as.numeric(x.fc[,"X14h_vs_ctrl"]))##z28 4h vs 0h 

# determine boundaries (not so important)
#eff.max.val <- min( abs( c( max(x.fc), min(x.fc) ) ) )
eff.max.val <- 5
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

names(fc.icos) <- x.fc[,1]
names(fc.bbz) <- x.fc[,1]
names(fc.z28) <- x.fc[,1]

fc <- cbind(fc.icos,fc.bbz,fc.z28)

#fl.nm <- "your file name here"
#pdf(fl.nm)
#m <- m.cut # here goes the fc data
for(j in 2:ncol(x)){
	nm <- colnames(x)[j]
	y <- x[order(x$Th1,decreasing = F),] ## take non-zero rows 
	gns <- y[,1][which(y[,j]>0)]
	
	# cluster the rows
	#d <- dist(as.matrix(fc[gns,]))
	d <- as.matrix(fc[gns,])
	temp <- d[,1]
	d[,1] <- d[,3]
	d[,3] <- temp ## change the order of the columns
	colnames(d) <- c("28z","BBz","ICOSz")
	#hc <- sort(d, decreasing = T)
	#gn.order <- hc$labels[hc$order]
	gn.order <- gns
	
# plot heatmap
	if(length(gns)>1){
	                pheatmap(d,scale="none",color=colorRampPalette(c(orange.col,"white",
	                blue.col))(255),legend_breaks = c(-5,0,5),breaks=seq(-1*eff.max.val,eff.max.val,by=(2*eff.max.val)/255),cellwidth = 36, cellheight = 12,treeheight_col=0,
	                cex=1,cluster_cols=F,cluster_rows=F,main=colnames(x)[j])

	        }
#	pheatmap(x.fc[gn.order],scale="none",cellheight = 4,main="FC",breaks=seq(-eff.max.val,eff.max.val,by=(2*eff.max.val)/tot ),
#		colorRampPalette(c(orange.col,"white",blue.col))(tot),	
#	cex=0.3,main=paste(sep="","kc=",cut.abs," pval=",deseq.pval.cut),cluster_rows=F,cluster_cols=F)
	
}

a <- c(rep(0,128),rep(1,128))
b <- c(seq(0,1,by=1/255))
for(i in 1:256) 
{ 
polygon(c(-0.05, -0.05, 0.05, 0.05), c(b[i], b[i+1], b[i+1], b[i]), border=0, col = colorRampPalette(c(orange.col,"white",blue.col))(256)[i]) 
} 
 
dev.off()