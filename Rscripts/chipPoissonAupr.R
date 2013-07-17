file <- '~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_cent_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,4]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.cent.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.cent.luc$AUPR
##########
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_htseq_format_cut' >20
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,4]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.cut.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.cut.luc$AUPR
########## >20
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_htseq_gname_format_cut'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,3]
pred.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.gname.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.gname.luc$AUPR
###########
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_YY1_hg19_htseq_poisson_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,2],pred.gs.yy1[,3],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.poisson.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.poisson.luc$AUPR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_cent_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,4]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.cent.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.cent.luc$AUPR
##########
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_htseq_format_cut' #>20
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,4]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- paste(pred.gs.yy1[,1],pred.gs.yy1[,2],sep= ":")
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.cut.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.cut.luc$AUPR
########## 
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_htseq_gname_format_cut'#>20
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,3]
pred.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.gata2.htseq.gname.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.gata2.htseq.gname.luc$AUPR
###########
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_htseq_poisson_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,2],pred.gs.yy1[,3],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.gata2.htseq.poisson.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.gata2.htseq.poisson.luc$AUPR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Gata2sc267Pcr1xPkRep1genes_poisson_prox_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.gata2.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.gata2.htseq.poisson.chip$AUPR
###################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Gata2sc267Pcr1xPkRep1genes_gname_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,2]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.gata2.htseq.gname.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.gata2.htseq.gname.chip$AUPR
####################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562MaxV0416102PkRep1genes_poisson_prox_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.max.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.max.htseq.poisson.chip$AUPR
###################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562MaxV0416102PkRep1genes_gname_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,2]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.max.htseq.gname.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.max.htseq.gname.chip$AUPR
####################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Yy1sc281V0416101PkRep1genes_poisson_prox_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.poisson.chip$AUPR
###################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Yy1sc281V0416101PkRep1genes_gname_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,2]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.gname.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.gname.chip$AUPR
####################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562CtcfcPcr1xPkRep2genes_poisson_dist_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.ctcf.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.ctcf.htseq.poisson.chip$AUPR
#####################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562NrsfV0416102PkRep1genes_poisson_dist_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.nrsf.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.nrsf.htseq.poisson.chip$AUPR
####################
plot(aupr.max.htseq.poisson.chip$rec, aupr.max.htseq.poisson.chip$prec, col = 'skyblue', type = 'lines', xlab = 'recall', ylab = 'precision',lwd=3,ylim=c(0,1))
lines(aupr.yy1.htseq.poisson.chip$rec, aupr.yy1.htseq.poisson.chip$prec, col = 'royalblue', type = 'lines',lwd=3)
lines(aupr.gata2.htseq.poisson.chip$rec, aupr.gata2.htseq.poisson.chip$prec, col = 'red', type = 'lines',lwd=3)
lines(aupr.nrsf.htseq.poisson.chip$rec, aupr.nrsf.htseq.poisson.chip$prec, col = 'deepskyblue', type = 'lines',lwd=3)
lines(aupr.ctcf.htseq.poisson.chip$rec, aupr.ctcf.htseq.poisson.chip$prec, col = 'orange', type = 'lines',lwd=3)
legend("topright",bty='n',lty=1,lwd=3,legend=c(paste(sep="","AUPR.MAX=",round(aupr.max.htseq.poisson.chip$AUPR,2) ), paste(sep="","AUPR.YY1=",round(aupr.yy1.htseq.poisson.chip$AUPR,2) ),paste(sep="","AUPR.NRSF=",round(aupr.nrsf.htseq.poisson.chip$AUPR,2) ),
paste(sep="","AUPR.CTCF=",round(aupr.ctcf.htseq.poisson.chip$AUPR,2) ),paste(sep="","AUPR.GATA2=",round(aupr.gata2.htseq.poisson.chip$AUPR,2) )),col=c("skyblue","royalblue","deepskyblue","orange","red"))
