if (!require(tidyverse, quietly = T)) install.packages("tidyverse")
if (!require(adegenet, quietly = T)) install.packages("adegenet")
if (!require(ape, quietly = T)) install.packages("ape")
if (!require(seqinr, quietly = T)) install.packages("seqinr")
if (!require(data.table, quietly = T)) install.packages("data.table")
if (!require(data.table, quietly = T)) install.packages("xlsx")

library(tidyverse)
library(adegenet)
require(ape)
require(seqinr)
require(data.table)
require(xlsx)
library (plyr)

#setwd("G:/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Database/phylogenetic signal")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Database/phylogenetic signal")

phylosygnal<-function(file,color,name){
  
  #Neighbour-Joining tree (nj) based on our new distance matrix.
  
  labels <- labels(file)
  #d_split = strsplit(labels, "[|]")
  #fasta_align_names = data.frame(t(sapply(d_split,`[`)))
  #fasta_align_names$labels <- labels
  #names(fasta_align_names)<-c("strain","date","region","country","state","labels")
  
  dates = as.Date(sapply(strsplit(labels(file), "[|]"), function(ele) ele[2]))
  
  fasta_align_names = data.frame(labels,dates)
  
  dates <- as.Date(fasta_align_names$date)
  head(dates)
  
  #fasta_align_names[31,]
  range <- range(dates)
  range
  
  #which(fasta_align_names$date == "2019-12-30", arr.ind = TRUE)
  
  days <- as.integer(difftime(dates, min(dates), unit="days"))
  days
  
  ### DNA distances
  
  D.tn93 <- dist.dna(file, model="TN93")
  
  tre <- nj(D.tn93)
  
  ## distances from the root
  D.from.root <- cophenetic(tre)[1,]
  
  ## plot distances vs time
  
  #Time from the root:days
  
  plot(days, D.from.root, pch=16, col=transp(color), cex=5, cex.lab=1.5, cex.axis=1.5,
       xlab="Days", ylab="Distance from the root", main=name, cex.main=5)
  
  
  ## calibrate clock (linear regression through origin)
  
  clock <- lm(D.from.root~-1+days)
  clock
  abline(clock, col=color)
  return(plot)
}





#####################

#### B_1

B_1 = fasta2DNAbin("B1_nextclade.aligned.fasta")
B_1 = B_1[, 56:29803]

#### B_1_111

B_1_111 = fasta2DNAbin("B1_111_nextclade.aligned.fasta")
B_1_111 = B_1_111[, 56:30181]

#### B_1_348

B_1_348 = fasta2DNAbin("aligned_B1_1_348_downsampling.fasta")
B_1_348 = B_1_348[, 56:29804]

#### B_1_420

B_1_420 = fasta2DNAbin("B.1.420_downsampling_nextclade.aligned.fasta")
B_1_420 = B_1_348[, 56:29803]

#### P_1
P_1 = fasta2DNAbin("P_1_nextclade.aligned.fasta")
P_1 = P_1[, 56:29803]

#### B_1_1_7

B_1_1_7 = fasta2DNAbin("aligned_B1_1_7_downsampling.fasta")
B_1_1_7 = B_1_1_7[, 56:29817]

#### C_37

C_37 = fasta2DNAbin("aligned_C_37_downsampling.fasta")
C_37 = C_37[, 56:29807]


#Delta 

Delta = fasta2DNAbin("delta_downsampling_nextclade_aligned.fasta")
Delta = Delta[, 56:29803]

#Omicron 

Omicron = fasta2DNAbin("Omicron_downsampling_nextclade_aligned.fasta")
Omicron = Omicron[, 56:29803]

pdf("S1Fig.pdf",width=20,height=10)
par(mfrow=c(3,3))

phylosygnal(B_1,"#91D1C2FF","B.1")
phylosygnal(B_1_111,"#64CC80","B.1.111")
phylosygnal(B_1_348,"#80A680","B.1.1.348")
phylosygnal(B_1_420,"magenta","B.1.420")
phylosygnal(C_37,"#3C5488FF","Lambda")
phylosygnal(B_1_1_7,"#DC0000FF","Alpha")
phylosygnal(P_1 ,"#00A087FF","Gamma")
phylosygnal(Delta ,"#4DBBD5B2","Delta")
phylosygnal(Omicron ,"#8491B4FF","Omicron")

dev.off()
