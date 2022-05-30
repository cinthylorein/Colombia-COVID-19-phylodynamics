library(bdskytools)
library(ggplot2)
library(pammtools)
library(ggpubr)
library(tidyverse)
library(adegenet)
library(ape)
library(data.table)
library(seqinr)
library("scales")


p_rt<-function(file1,file2,file3,file4,file5,file6,file7,file8,file9,file10){
  B_1<-as.data.frame(fread(file1))
  B_1_111<-as.data.frame(fread(file2))
  B_1_1_348<-as.data.frame(fread(file3))
  C37<-as.data.frame(fread(file4))
  B_1_1_7<-as.data.frame(fread(file5))
  P1<-as.data.frame(fread(file6))
  B_1_621<-as.data.frame(fread(file7))
  Delta<-as.data.frame(fread(file8))
  Omicron<-as.data.frame(fread(file9))
  B_1_420<-as.data.frame(fread(file10))
  #min(B$date)
  #setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")
  #setwd("G:/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")
  #Mob <- as.data.frame(fread("Col_mobility.csv"))
  #Mob$mobility_mean <- abs(Mob$mobility_mean)
  #B <- merge(B,Mob, by.x="date", by.y="date")
  #cases <- as.data.frame(fread("descriptive_database_OWID.tsv.tsv"))
  #cases$new_cases <- log(cases$new_cases)
  #B <- merge(B,cases, by.x="date", by.y="date")
  #name<-sapply(strsplit(file1,""),"[",1)
  coeff <- 1
  plot<-ggplot(B_1,aes(x=time,y=mean))+
    theme_set(theme_classic(base_size=20))+
    theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), 
          axis.title = element_text(size = 20), axis.title.y = element_text(face = "italic"), 
          plot.title = element_text(face="bold", size=20))+
    geom_step(color="#91D1C2FF")+
    #geom_stepribbon(aes(ymin=upper,ymax=lower),alpha=0.9,fill=color)+
    geom_line(data= B_1_111, aes(x=time,y=mean), color="#64CC80")+
    geom_line(data= B_1_1_348, aes(x=time,y=mean), color="#80A680")+
    geom_line(data= C37, aes(x=time,y=mean), color="#3C5488FF")+
    geom_line(data= B_1_1_7, aes(x=time,y=mean), color="#DC0000FF")+
    geom_line(data= P1, aes(x=time,y=mean), color="#00A087FF")+
    geom_line(data= B_1_621, aes(x=time,y=mean), color="#F39B7FFF")+
    geom_line(data= Delta, aes(x=time,y=mean), color="#4DBBD5B2")+
    geom_line(data= Omicron, aes(x=time,y=mean), color="#8491B4FF")+
    geom_line(data= B_1_420, aes(x=time,y=mean), color="magenta")+
    #geom_stepribbon(data= B, aes(ymin=upper,ymax=lower), alpha=0.2, color="grey")+
    #geom_line(data= B, aes(x=time,y=new_cases_p), lty=2)+
    #scale_y_continuous( name = "Ne", sec.axis = sec_axis(trans = ~.*coeff, name="% New cases"))+
    geom_vline(xintercept =2020.23, colour = "blue", lty=2)+ 
    # Colombia first case 2020-02-26 (31+26)/365=0.15
    geom_vline(xintercept =2020.15, colour = "black", lty=1)+
    #wave1 Aug 2020 214/365
    geom_vline(xintercept =2020.58, colour = "grey", lty=1, size=3, alpha=0.5)+
    #wave2 Jan 2021  25/365
    geom_vline(xintercept =2021.06, colour = "grey", lty=1, size=3, alpha=0.5)+
    #wave3 Jun 2021 152/365
    geom_vline(xintercept =2021.41, colour = "grey", lty=1, size=3, alpha=0.5)+
    #wave4 Jan 2022
    geom_vline(xintercept =2022.005, colour = "grey", lty=1, size=3, alpha=0.5)+
    labs(x="Date", y="Ne")+
    xlim(c(2020,2022.5))
    #ylim(c(0,150))+
    #ggtitle(name)
    return(plot)
}

#theme_set(theme_bw(base_size=12))
setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")
#setwd("G:/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")


plot <- p_rt("COVID_BICEPS_B_1","COVID_BICEPS_B_1_111","COVID_BICEPS_B_1_1_348","COVID_BICEPS_C_37","COVID_BICEPS_B_1_1_7","COVID_BICEPS_P_1","COVID_BICEPS_B_1_621","COVID_BICEPS_Delta","COVID_BICEPS_Omicron","COVID_BICEPS_B_1_420")


setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure3")
pdf ('Figure3_D.pdf',width=10,height=10)
#par(mfrow=c(2,4))
plot 
dev.off()
