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
aupr.yy1.htseq.gname.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.gname.luc$AUPR
###########
file <- '~/Documents/google-python-exercises/gold_standard/luciferase_GATA2_hg19_htseq_poisson_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,2],pred.gs.yy1[,3],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,4]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.poisson.luc <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.poisson.luc$AUPR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Gata2sc267Pcr1xPkRep1genes_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.poisson.chip$AUPR
###################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Gata2sc267Pcr1xPkRep1genes_gname_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,2]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.gname.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.gname.chip$AUPR
####################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562MaxV0416102PkRep1genes_format'
pred.gs.yy1 <- read.delim(file,header=F)
pred.yy1 <- -log10(ppois(pred.gs.yy1[,3],pred.gs.yy1[,4],lower.tail=F))
gs.yy1 <- pred.gs.yy1[,2]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.poisson.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.poisson.chip$AUPR
###################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562MaxV0416102PkRep1genes_gname_format'
pred.gs.yy1 <- read.delim(file,header=F)
gs.yy1 <- pred.gs.yy1[,2]
pred.yy1 <- pred.gs.yy1[,3]
names(gs.yy1) <- pred.gs.yy1[,1]
names(pred.yy1) <- names(gs.yy1)
gs.yy1 <- gs.yy1>0
aupr.yy1.htseq.gname.chip <- calcAupr(pred.yy1,gs.yy1)
aupr.yy1.htseq.gname.chip$AUPR
####################
file <- '~/Documents/google-python-exercises/gold_standard/wgEncodeHaibTfbsK562Yy1sc281V0416101PkRep1genes_format'
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
plot(aupr.yy1.htseq.poisson.chip$rec, aupr.yy1.htseq.poisson.chip$prec, col = 'slategray', type = 'lines', xlab = 'recall', ylab = 'precision')
lines(aupr.yy1.htseq.poisson.chip$rec, aupr.yy1.htseq.poisson.chip$prec, col = 'royalblue', type = 'lines')
lines(aupr.yy1.htseq.poisson.chip$rec, aupr.yy1.htseq.poisson.chip$prec, col = 'skyblue', type = 'lines')
legend("topright",lty=1,lwd=2,legend=c(paste(sep="","AUPR.MAX=",0.78 ), paste(sep="","AUPR.YY1=",0.73 ),paste(sep="","AUPR.GATA2=", 0.20)),col=c("royalblue","slategray","skyblue"))
