if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
install.packages("ggplot2")
install.packages("ggimage")
install.packages("igraph")
install.packages("treeio")

library(ggtree)
library(ggplot2)
library(ggimage)
library("igraph")
library(treeio)
library(ggpubr)


tree<-function(file1,name, se1){
  beast <- read.beast(file1)
  ggtree(beast, mrsd=se1, aes(color=location)) +theme_tree2()+  
    #geom_nodelab(aes(x=branch, label=round(as.numeric(posterior), 2)), vjust=-.5, size=2) +
    geom_range(range='height_0.95_HPD', color='grey', alpha=.6, size=1) +
    geom_tippoint(aes(color=location), size=1, alpha=.75) +
    geom_nodepoint(aes(size=as.numeric(posterior)), alpha=.1)+
  #  scale_color_brewer("location", palette="viridis") +
    geom_hilight(node=50, fill="gold")+
    theme_tree2(legend.position = "right")+ 
    scale_y_continuous(expand=c(0, 0.6)) + 
#    layout_dendrogram()+
    ylab(name)
}

setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Colombia_RW/BICEPS_DTA/B_1")
B_1<-tree("COVID_DTA_ANN.tree","B.1","2021-07-21")
setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Colombia_RW/BICEPS_DTA/B_1_111")
B_1_111<-tree("COVID_DTA_ANN.tree","B.1.111","2021-08-27")
setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Colombia_RW/BICEPS_DTA/B_1_1_348")
B_1_1_348<-tree("COVID_DTA_ANN.tree","B.1.1.348","2021-05-18")
setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Colombia_RW/BICEPS_DTA/P_1")
P1 <-tree("COVID_DTA_ANN.tree","Gamma","2021-08-25")
setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Colombia_RW/BICEPS_DTA/B_1_1_7")
B_1_1_7 <- tree("COVID_DTA_ANN.tree","B.1.1.7","2021-05-25")
setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Colombia_RW/BICEPS_DTA/C_37")
C37 <- tree("COVID_DTA_ANN.tree","Lambda","2021-08-23")

Delta <- tree("COVID_DTA_ANN.tree","Lambda","2021-08-23")

p_all<-ggarrange(B_1,P1,B_1_111,B_1_621,B_1_1_348,B_1_1_7,C37, ncol = 2, nrow = 4)

setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Figures")
ggsave("S3Fig.jpg",plot=p_all,dpi=500,width=12,height=15)


