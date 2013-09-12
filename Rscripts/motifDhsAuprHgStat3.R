rm(list=ls())
## gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"
# load("/Users/xichen/data/results_combine_data_exprsn_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1.Rdata")
# #pred.file <- 
# #knockout_rna_zscores.all <- read.delim(gs.file)
# 
# stat3.gs <- res[["activation"]][["IRF4"]][["th17"]][,1:15]
# #stat3_gs <- cbind(stat3_kcri,stat3_kc,stat3_kri,stat3_k,stat3_cri,stat3_c,stat3_ri)
# #+res[["activation"]][["BATF"]][["th17"]][,"KCRI"]+res[["activation"]][["IRF4"]][["th17"]][,"KCRI"]+res[["activation"]][["MAF"]][["th17"]][,"KCRI"]+res[["activation"]][["RORC"]][["th17"]][,"KCRI"]
# 
# mm.to.hg <- read.delim("~/data/HMD_HumanSequence.txt", header=F)
# mm.to.hg <- as.matrix(mm.to.hg)
# 
# stat3.gs.hg <- matrix(NA,nr = length(intersect(rownames(stat3.gs),toupper(mm.to.hg[,1]))),nc = ncol(stat3.gs))
# k <- 1
# rownames(stat3.gs.hg) <- rep(NA,length(intersect(rownames(stat3.gs),toupper(mm.to.hg[,1]))))
# colnames(stat3.gs.hg) <- colnames(stat3.gs)
# for (i in 1:nrow(stat3.gs)){
#     for (j in 1:nrow(mm.to.hg)){
# 	if (toupper(mm.to.hg[j,1]) == rownames(stat3.gs)[i]){
# 	    stat3.gs.hg[k,] <- stat3.gs[i,]
# 	    rownames(stat3.gs.hg)[k] <- toupper(mm.to.hg[j,4])
#             k <- k+1
# 	}
#     }
# }
# 
# save.image('~/code/AUPR/th17Stat3MmToHg.RData')

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
##AUPR script to compare prediction to gold standard of different sizes
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

source('~/code/AUPR/aupr.R')
debug <- TRUE

print.error <- function(){
	cat("
DESCRIPTIION:	
	write something?

INPUT:
	1.input_files: path to MACS/bed files '::' delim [path_input=f1::f2::f3::...::fk]
	2.path_output: path to save generated MTL cluster file (where to save mtls.xls)
	3.expt_names: user specified names for MACS files '::' delim [expt_names=n1::n2::n3::...::nk]
	4.input_type: the type of input file used (MACS or BED; defaults to MACS)
	5.mtl_type: interval or summit (defaults to summit)
	6.dist.summits: maximum distance between summits belonging to the same MTL (defaults to 100; only used if mtl_type is summit)
		
EXAMPLE RUN: 
	cluster_peaks.R
	--input_files input/SL2870_SL2871_peaks.xls::input/SL2872_SL2876_peaks.xls::input/SL3032_SL2871_peaks.xls::input/SL3037_SL3036_peaks.xls::input/SL3315_SL3319_peaks.xls
	--input_type MACS
	--path_output results/
	--expt_names RORC_Th17::IRF4_Th17::MAF_Th17::BATF_Th17::STAT3_Th17
	--dist_summits 100
	--mtl_type summit\n\n")
}

### retrieve command line args or enter into debug mode
if(debug==T){
	cmd.args <- c(
		#"--input_gold_standard_files ~/code/AUPR/th17Stat3MmToHg.RData"#"data/xls/SL10571_SL10565_peaks.xls::data/xls/SL10570_SL10564_peaks.xls::data/xls/SL10572_SL10566_peaks.xls",
		"--input_prediction_files ~/data/zscore/STAT3_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/IRF4_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/BATF_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/RORC_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/FOSL2_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/GATA3_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/TBX21_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/CTCF_DGF_Th17_htseq_poissoncdf_zscore::~/data/zscore/RANDOM_DGF_Th17_htseq_poissoncdf_zscore",
		#"--input_type MACS", BED
		"--path_output ~/data/AUPR/STAT3/",
		"--tf_name STAT3::IRF4::BATF::RORC::FOSL2::GATA3::TBX21::CTCF::RANDOM",
		"--gs_name STAT3"
		#"--z_cutoff 0"
		#"--expt_names macs1::macs2::macs3",
		#"--expt_names bed1::bed2::bed3",
		#"--mtl_type interval", #interval summit
		#"--dist_summits 100"
	)
} else {
	cmd.args <- commandArgs(trailingOnly = T);      
}

args.nms.must <- c(	
		#"--input_gold_standard_files",   #1
		"--input_prediction_files",
		"--path_output",  	 #2
		"--tf_name",
		"--gs_name"
)

# read command line paramters that are not optional
read.cmd.line.params.must <- function(args.nms, cmd.line.args){
	#if(length(grep("--version",cmd.line.args))){
	#	cat("version",script.version,"\n")
	#	q()
	#}
	args <- sapply(strsplit(cmd.line.args," "),function(i) i)
	vals <- character(length(args.nms))
	# split cmd.line to key and value pairs
	for(i in 1:length(args.nms)){
		ix <- grep(args.nms[i],args)
		if(length(ix)>1){
			stop("arg ",args.nms[i]," used more than once.  Bailing out...\n",print.error())
		} else if (length(ix)==0){
			stop("could not find ",args.nms[i],". Bailing out...\n",print.error())
		} else {
			vals[i] <- args[ix+1]
		}
	}
	return(vals)
}

vals.must <- read.cmd.line.params.must(args.nms = args.nms.must, cmd.line.args = cmd.args)
#input.gs.files <- vals.must[1]
input.pred.files <- vals.must[1]
input.pred.files.list <- strsplit(input.pred.files,"::")[[1]]
path.output <- vals.must[2]
tfnames <- vals.must[3]
tfnames.list <- strsplit(tfnames,"::")[[1]]
gsname <- vals.must[4]

nlist <- seq(50,1000,by=50)

pdf(file=paste(sep="",path.output,"aupr_",gsname,".pdf"))
plot(nlist,rep(0,length(nlist)),lwd = 1,xlab = paste(sep="","Length of ",gsname,"-net gold standard"),ylab = 'AUPR',type ='l',ylim=c(0,0.1),col='black',xaxt="n")
axis(1, at = nlist, las=2)
########################################################################################################################

load('~/code/AUPR/th17Stat3MmToHg.RData')
gs.hg <- stat3.gs.hg
gs.hg.df <- as.data.frame(gs.hg)
gs.hg.kcri <- gs.hg.df[order(gs.hg.df$C,decreasing=TRUE),]
#expt.names <- strsplit(vals.must[3],"::")[[1]]

#if(length(input.files)==1){
#	cat("only provided one MACS file to cluster.")
#	print.error()
#}

col.list <- c('navy','royalblue','skyblue','lightblue','paleturquoise','red','brown','orange','gray')
## loop over the number of files
for (f in 1:length(tfnames.list)){
	input.pred.file <- input.pred.files.list[f]
	tfname <- tfnames.list[f]
	#input.pred.files <- '~/hg19priors/v081913/STAT3_DGF_Th17_htseq_poissoncdf_zscore'
	motif.and.dnase.reads.stat3.rep1 <- read.delim(input.pred.file, header=F)
	#motif.and.dnase.reads.stat3.rep1 <- read.delim("~/data/zscore/IRF4_DGF_Th17_htseq_output_poissoncdf_zscore", header=F)

	m_d_network <- motif.and.dnase.reads.stat3.rep1[,2]
	#hist(m_d_network_stat3)
	names(m_d_network) <- toupper(motif.and.dnase.reads.stat3.rep1[,1])

	aupr_res <- matrix(NA, nr = length(nlist), nc = 3)
	auroc_res <- matrix(NA, nr = length(nlist), nc = 1)

	for (n in 1:length(nlist)){
	        gs.hg.kcri.cut <- gs.hg.kcri[1:nlist[n],]
	        unique.names <- unique(c( names(m_d_network), rownames(gs.hg.kcri.cut)))
	        overlap <- (length(m_d_network) + length(gs.hg.kcri.cut[,'C'])) - length(unique.names)
	        data <- matrix(NA, nr=length(unique.names), nc=2)
	        rownames(data) <- unique.names
	        colnames(data) <- c("motif_dnase",'C')
	        for (i in 1:nrow(data)) {
	        data[i,1] <- m_d_network[unique.names[i]]##read count per gene
	        if (unique.names[i] %in% rownames(gs.hg.kcri.cut)){
	            data[i,2] <- gs.hg.kcri.cut[unique.names[i],'C']##kcri network score
	        }
	    }
	    pred <- data[,1]
	    pred <- pred[which(!is.na(pred))]
	    gs <- data[,2]
	    gs[which(is.na(gs))] <- 0
	    gs <- gs>0
	    aupr <- calcAupr(pred,gs)
	    aupr_res[n,1] <- aupr$AUPR
	        auroc_res[n,1] <- aupr$AUROC
	        aupr_res[n,2] <- overlap
	        aupr_res[n,3] <- length(unique.names)
	}
	rownames(aupr_res) <- nlist
	colnames(aupr_res) <- c('C','intersect','union')
	rownames(auroc_res) <- nlist
	colnames(auroc_res) <- c('C')

	write.table(aupr_res,sep="\t",,file=paste(sep="",path.output,"aupr_",tfname,"_",gsname))
	write.table(auroc_res,sep="\t",,file=paste(sep="",path.output,"auroc_",tfname,"_",gsname))
	
	lines(nlist,aupr_res[,'C'],col=col.list[f],lwd=2)
	
}

legend("topleft",lty=1,lwd=2,legend=tfnames.list,col=col.list)
dev.off()



#icosp_aupr <- calcAupr(pred,gs)
#plot(icosp_aupr$rec, icosp_aupr$prec, type = 'lines', main = 'IRF4_MOTIF1E3_DGF_VS_KCRI', xlab = 'recall', ylab = 'precision')
#lines(bb_aupr$rec, bb_aupr$prec, col="darkgreen")
#lines(z28_aupr$rec, z28_aupr$prec, col="blue")
#legend("topright",lty=1,lwd=2,legend=c("AUPR=0.6955647","AUROC=0.4749878"))

#save.image('~/code/auprVgs.RData')
