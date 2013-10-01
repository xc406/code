##python ~/code/htseq/faireMat.py Sample_lane5.E14.sorted.bam Sample_lane6.E14.sorted.bam Sample_lane1.EpiSC.sorted.bam Sample_lane3.EpiSC.sorted.bam Sample_lane1.NSC.sorted.bam Sample_lane2.NSC.sorted.bam Sample_lane7.MEF.sorted.bam Sample_lane8.MEF.sorted.bam Sample_lane7.gDNA.sorted.bam Sample_lane8.gDNA.sorted.bam genome1kb.gff 1000 27160662 27283392 30836626 28698142 30161179 32705690 28374688 28330417 25950674 25326527

d <- read.delim('~/Documents/google-python-exercises/fireMatrixECleanExpr',header=T)
d.df <- d[,2:13]#26]#652]
rownames(d.df) <- d[,1]
for (i in length(colnames(d.df))){
	if (strsplit(colnames(d.df)[i],"[.]")[[1]][1] %in% mylist){
		mylist <- mylist
	}else{
	    mylist[length(mylist)+1] <- strsplit(colnames(d.df)[i],"[.]")[[1]][1]
	}
}
for (i in colnames(d.df)){
	for (j in colnames(d.df)){
		if (strsplit(colnames(d.df)[i],"[.]")[[1]][1]==strsplit(colnames(d.df)[j],"[.]")[[1]][1]){
			d.df[,i] = d.df[,i]+d.df[,j]
			drops <- c(colnames(d.df)[j])
			d.df <- d.df[,!(colnames(d.df) %in% drops)]
		}
	}
	
}
d.df[,'Sp4'] = d.df[,'Sp4']+d.df[,'Sp4.1']
d.df[,'Egr1'] = d.df[,'Egr1']+d.df[,'Egr1.1']
d.df[,'Pou5f1'] = d.df[,'Pou5f1']+d.df[,'Pou5f1.1']+d.df[,'Pou5f1.2']+d.df[,'Pou5f1.3']
d.df[,'Sox2'] = d.df[,'Sox2']+d.df[,'Sox2.1']+d.df[,'Sox2.2']
d.df[,'Klf4'] = d.df[,'Klf4']+d.df[,'Klf4.1']
d.df[,'Zfp281'] = d.df[,'Zfp281']+d.df[,'Zfp281.1']
drops <- c("Sp4.1","Egr1.1","Pou5f1.1","Pou5f1.2","Pou5f1.3","Sox2.1","Sox2.2","Klf4.1","Zfp281.1")
d.df.clean <- d.df[,!(colnames(d.df) %in% drops)]
d.scaled.clean <- as.matrix(scale(d.df.clean))

d.scaled <- as.matrix(scale(d.df))
drops <- c("M1794_Prrx2","M1758_Gata2","M1525_Tbp..Tbpl2","M0882_Lhx4","M0939_Lmx1b","M0956_Prop1","M1807_Sox5","M1784_Nkx2.5","M1863_Arid3a","M1866_Nfic","M1908_Six5","M1845_Pdx1")
drops <- c("M1794_Prrx2","M1758_Gata2")
d.scaled.clean <- d.scaled[,!(colnames(d.scaled) %in% drops)]
##heatmap(d.scaled.clean, scale='none')##too many elements

hc.fire <- hclust(dist(d.scaled))
hc.motif <- hclust(dist(t(d.scaled)))
#heatmap(d.scaled.clean[cutree(hc.fire,k=5)==1,],Colv=as.dendrogram(hc.motif), scale='none')
heatmap(d.scaled.clean[cutree(hc.fire,k=1000)==1,cutree(hc.motif,k=30)==1], scale='none')
library("RColorBrewer")
color.ramp<-colorRampPalette(c("black","darkblue", "blue","yellow","orange","red"))(50)
heatmap(d.scaled.clean[cutree(hc.fire,k=5)==1,], scale='none',labRow=F,cexCol = 1,col = color.ramp)
                                 
or.col <- c(255,102,0)
bl.col <- c(0,51,255)
blue.col <- rgb(bl.col[1],bl.col[2],bl.col[3],max=255)
orange.col <- rgb(or.col[1],or.col[2],or.col[3],max=255)      
color.ramp <- colorRampPalette(c(orange.col,rep("white",4),blue.col))(255)                            
heatmap(d.scaled.clean[cutree(hc.fire,k=4)==1,], scale='row',labRow=F,cexCol = 1,Colv=FALSE,Rowv=FALSE)     
##fire matrix heatmap

color.ramp <- colorRampPalette(c("ivory","lemonchiffon","red"))(255)
heatmap(d.scaled, scale=none,labRow=F,cexCol = 1,Colv=FALSE,Rowv=FALSE,col=color.ramp)
##Jaccard Index heatmap
d.matrix <- as.matrix(d.df)
heatmap(d.matrix, Rowv = FALSE,Colv = FALSE,revC=TRUE,col=color.ramp)
heatmap(d.matrix, Rowv = FALSE,Colv = FALSE,revC=TRUE,scale='none')
meme esc_fire_dist.fa -dna -mod anr -revcomp -bfile mm9bgfile  


nP <- list(col = 3:2, pch =  21:22, lab.cex = 0.5, lab.col = "darkgray")
plot(as.dendrogram(hc.motif),horiz=TRUE,edgePar = list(col = "gray", lwd = 1),nodePar= nP)      
heatmap(d.scaled[cutree(hc.fire,k=4)==1,], scale='row',cexCol = 1,labRow=F,Colv=FALSE,Rowv=FALSE)
heatmap(d.scaled, scale='row',cexCol = 1,labRow=F,Colv=FALSE,Rowv=FALSE)

fire.order <- hc.fire$labels[hc.fire$order] ##return fire elements in the clustered order
write.table(fire.order,sep='\n',quote=FALSE,row.names=FALSE,col.names=FALSE,file='~/Documents/google-python-exercises/fireOrder')

#fire.order <- read.delim('~/Documents/google-python-exercises/fireOrder',header=F)
#oct4.plot <- data.frame(matrix(NA, nr=length(fire.order), nc=9))
for (i in 1:length(fire.order)){
	fire.order[i] <- paste(strsplit(fire.order[i],'_')[[1]][1],strsplit(fire.order[i],'_')[[1]][2],strsplit(fire.order[i],'_')[[1]][3],sep='_')
#	for (j in 1:nrow(oct4.df)){
#		if (length(grep(id[i], rownames(oct4.df)[j]))){
##			#m <- 5
#			oct4.plot[i,] <- oct4.df[j,]
##			#rownames(oct4.plot[i,]) <- id 
#		}
#	}
}

oct4 <- read.delim('~/Documents/google-python-exercises/chipplotSox2e',header=T)
oct4.df <- oct4[,1:100]
#rownames(oct4.df) <- oct4[,1]#paste(strsplit(oct4[,1],'_')[[1]][1],strsplit(oct4[,1],'_')[[1]][2],strsplit(oct4[,1],'_')[[1]][3],sep='_')#
colnames(oct4.df) <- seq(1,100,by=1)
for (i in 1:nrow(oct4.df)){
	rownames(oct4.df)[i] <- paste(strsplit(rownames(oct4.df)[i],'_')[[1]][1],as.character(as.numeric(strsplit(rownames(oct4.df)[i],'_')[[1]][2])+1),strsplit(rownames(oct4.df)[i],'_')[[1]][3],sep='_')
	#rownames(oct4.df)[i] <- paste(strsplit(rownames(oct4.df)[i],'_')[[1]][1],strsplit(rownames(oct4.df)[i],'_')[[1]][2],strsplit(rownames(oct4.df)[i],'_')[[1]][3],sep='_')
}

oct4.df.sorted <- oct4.df[match(fire.order,rownames(oct4.df)),]##order df by fire.order
oct4.df.sorted[oct4.df.sorted==0.0] <- 0.000001
oct4.df.sorted.log2 <- log2(oct4.df.sorted)

wce <- read.delim('~/Documents/google-python-exercises/chipplotWCEe',header=T)
wce.df <- wce[,1:100]
#rownames(wce.df) <- wce[,1]#paste(strsplit(oct4[,1],'_')[[1]][1],strsplit(oct4[,1],'_')[[1]][2],strsplit(oct4[,1],'_')[[1]][3],sep='_')#
colnames(wce.df) <- seq(1,100,by=1)
for (i in 1:nrow(wce.df)){
	rownames(wce.df)[i] <- paste(strsplit(rownames(wce.df)[i],'_')[[1]][1],as.character(as.numeric(strsplit(rownames(wce.df)[i],'_')[[1]][2])+1),strsplit(rownames(wce.df)[i],'_')[[1]][3],sep='_')
	#rownames(oct4.df)[i] <- paste(strsplit(rownames(oct4.df)[i],'_')[[1]][1],strsplit(rownames(oct4.df)[i],'_')[[1]][2],strsplit(rownames(oct4.df)[i],'_')[[1]][3],sep='_')
}

wce.df.sorted <- wce.df[match(fire.order,rownames(wce.df)),]
wce.df.sorted[wce.df.sorted==0.0] <- 0.000001
wce.df.sorted.log2 <- log2(wce.df.sorted)

oct4wce.df.sorted<- log2(oct4.df.sorted/wce.df.sorted)

pheatmap(oct4.df.sorted.log,scale="none",color=colorRampPalette(c("ivory","lemonchiffon",
                     "red"))(255),legend_breaks = c(-10,0,10),cellwidth = 1, cellheight = 0.15,treeheight_col=0,
                 cex=1,cluster_cols=F,cluster_rows=F,main='Nanog',show_rownames=FALSE,show_colnames=FALSE,legend=FALSE) 


################cell type clustering#########################
d <- read.delim('~/Documents/google-python-exercises/faireMatchr141000.txt',header=F)
colnames(d) <- c('E14.lane5','E14.lane6','EpiSc.lane1','EpiSc.lane3','NSC.lane1','NSC.lane2','MEF.lane7','MEF.lane8','gDNA.lane7','gDNA.lane8')
d.scaled <- as.matrix(scale(d))
hc.ct <- hclust(dist(t(d.scaled)))
nP <- list(col = 3:2, pch =  21:22, lab.cex = 1, lab.col = "darkgray")
par(pin=c(5,3))
plot(as.dendrogram(hc.ct),horiz=TRUE,edgePar = list(col = "gray", lwd = 2),nodePar= nP)

#d.merge <- matrix(NA, nc = 5, nr = nrow(d))
#d.matrix <- as.matrix(d)
gDNA <- (d$gDNA.lane7 + d$gDNA.lane8)/2
E14 <-(d$E14.lane5 + d$E14.lane6)/2
EpiSc <- (d$EpiSc.lane1 + d$EpiSc.lane3)/2
NSC <-(d$NSC.lane1 + d$NSC.lane2)/2
MEF <- (d$MEF.lane7 + d$MEF.lane8)/2
d.merge <- data.frame(E14,EpiSc,NSC,MEF,gDNA)
d.merge.scaled <- as.matrix(scale(d.merge))
hc.merge.ct <- hclust(dist(t(d.merge.scaled)))
par(pin=c(5,3))
plot(as.dendrogram(hc.merge.ct),horiz=TRUE,edgePar = list(col = "gray", lwd = 2),nodePar= nP)

MEF.minus <- MEF-gDNA
NSC.minus <- NSC-gDNA
EpiSc.minus <- EpiSc-gDNA
E14.minus <- E14-gDNA
d.merge <- data.frame(E14.minus,EpiSc.minus,NSC.minus,MEF.minus)

##exclude region in chr14


