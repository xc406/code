##load in data
motif.and.dnase.stat3.test <- read.delim("~/data/intersect_test7_two", header=F)

p1 <- 10**(-1*(motif.and.dnase.stat3.test[,2]))##calculate pvals from log transformed pvals
names(p1) <- toupper(motif.and.dnase.stat3.test[,1])
p2 <- 10**(-1*(motif.and.dnase.stat3.test[,3]))
names(p2) <- toupper(motif.and.dnase.stat3.test[,1])
n.names <- names(p1)
source('~/code/combine.p.R')
np <- list()
#p.combined<-list()
for (i in 1:length(p1)){
  np[i] <- combine.test(p = c(p1[i],p2[i]), method = "z.transform")##combine motif p-vals with dhs peak calling p-vals
  #p.combined[i] <- (-1)*log10(as.numeric(np[i]))
}
  
#m_d_network_stat3 <- as.numeric(np)
#names(m_d_network_stat3) <- names(p1)
#pname <- "NONAME"

####################
##scripts to combine p vals from multiple motifs/dhs combined p vals for one gene
####################
j <- 1
k <- 1
nr <- length(unique(c(names(m_d_network_stat3))))
v <- matrix(NA, nrow = nr, ncol = length(unique(c(names(m_d_network_stat3)))))
#v <- matrix(NA)
unique.names <- unique(c(names(m_d_network_stat3)))
rownames(v) <- unique.names
v[k,j] <- m_d_network_stat3[1]
#names(v[k,j]) <- names(m_d_network_stat3[1])
pname <- names(m_d_network_stat3[1])
for (i in 2:length(m_d_network_stat3)) {
  if (k < nr) {
    if (names(m_d_network_stat3[i]) == pname) {
      #names(v[k]) <- names(m_d_network_stat3[i])
      #v[k] <- (v[k])*(m_d_network_stat3[i])
      #v[k,] <- c(v[k,],rep(NA,1))
      v[k,j+1] <- m_d_network_stat3[i]
      #names(v[k,j+1]) <- names(m_d_network_stat3[i])
      j <- j+1
      } else {
      j <- 1
      k <- k+1
      pname <- names(m_d_network_stat3[i])
      v[k,j] <- m_d_network_stat3[i]
      #names(v[k,j]) <- names(m_d_network_stat3[i])
      }
    }
  }  
p.combined <- rep(NA, nr)
for (i in 1:nr){
   nc <- which(!is.na(v[i,]))
   pc <- rep(NA, length(nc))
   for (j in 1:length(nc)){     
     pc[j] <- v[i,j]
   } 
  p.combined[i] <- combine.test(p = pc, method = "z.transform")
}
names(p.combined) <- names(v[,1])

#######################
##use the combined p vals to calculate AUPR
#######################
m_d_network_stat3 <- p.combined
names(m_d_network_stat3) <- names(p.combined)

gs.file <- "~/Desktop/score_combine_KC_matrix_zcut_2.5_k_r_i_activator_Nov_28_2011.xls"
knockout_rna_zscores.all <- read.delim(gs.file)
k_r_network_stat3 <- knockout_rna_zscores.all[,"STAT3"]
names(k_r_network_stat3) <- toupper(rownames(knockout_rna_zscores.all))
k_r_network_stat3 <- k_r_network_stat3[which(k_r_network_stat3>0)]

#unique.names <- unique(c( names(m_d_network_stat3), names(k_r_network_stat3)))
unique.names <- unique(c(names(m_d_network_stat3), names(k_r_network_stat3)))

overlap <- (length(m_d_network_stat3) + length(k_r_network_stat3)) - length(unique.names)

data <- matrix(NA, nr=length(unique.names), nc=2)
rownames(data) <- unique.names
colnames(data) <- c("motif_dnase", "knockout_rna")
for (i in 1:nrow(data)) {
  data[i,1] <- m_d_network_stat3[unique.names[i]]
  data[i,2] <- k_r_network_stat3[unique.names[i]]
}
data[which(is.na(data))] <- 0

source('~/code/AUPR/aupr.R')
gs <- data[,2]
colnames(data)
pred <- data[,1]
gs>0
gs <- gs>0
gs
calcAupr(pred,gs)
stat3_aupr <- calcAupr(pred,gs)
plot(stat3_aupr$rec, stat3_aupr$prec, type = 'lines', main = 'Stat3_ztransform')
#x <- data.frame(n = n.names, "p-motif" = p1, "p-dhs" = p2, "p-combined" = np)
#rownames(x) <- n.names
#colnames(x) <- c("p-motif", "p-dhs", "p-combined")
#for (i in 1:length(p1)) {
#data[i,1] <- p1[motif.and.dnase.stat3.test[i]]
#   x[i,1] <- p1[i]
#   x[i,2] <- p2[i]
#   x[i,3] <- np[i]
#   #data[i,2] <- p2[motif.and.dnase.stat3.test[i,1]]
# }
