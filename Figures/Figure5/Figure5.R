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
library(dataRetrieval)
library(dplyr) # for `rename` & `select`
library(tidyr) # for `gather`
library(ggsci)
require(xlsx)
library(ggpubr)


p_rt<-function(file1,file2,color,se1,se2,name){
  SK<-as.data.frame(fread(file1))
  B<-as.data.frame(fread(file2))
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
  plot<-ggplot(SK,aes(x=time,y=mean))+
    theme_set(theme_classic(base_size=18))+
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
          axis.title = element_text(size = 18), axis.title.y = element_text(face = "italic"), 
          plot.title = element_text(face="bold", size=18))+
    geom_step(color=color)+
    geom_stepribbon(aes(ymin=upper,ymax=lower),alpha=0.7,fill=color)+
    geom_line(data= B, aes(x=time,y=mean))+
    geom_stepribbon(data= B, aes(ymin=upper,ymax=lower), alpha=0.2, color="grey")+
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
    # Colombia first squence reported
    geom_vline(xintercept =se1, colour = "red", lty=1)+ #se1
    # Lineage first squence reported 
    geom_vline(xintercept =se2, colour = "red", lty=2)+ #se2
    #text(mean(x1),(b*max(x1)/2+a),"text",srt=0.2,pos=3)
    labs(x="Date", y="Ne")+
    xlim(c(2020,2022.1))+
    #ylim(c(0,150))+
    ggtitle(name)
    return(plot)
}



#theme_set(theme_bw(base_size=12))
#setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")
#setwd("G:/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")

B_1<-p_rt("B.1","COVID_BICEPS_B_1","#91D1C2FF",2020.19,2020.021,"B.1")
B_1_111<-p_rt("B.1.111","COVID_BICEPS_B_1_111","#64CC80",2020.2,2020.18,"B.1.111")
B_1_1_348<-p_rt("B.1.1.348","COVID_BICEPS_B_1_1_348","#80A680",2020.33,2020.33,"B.1.1.348")
B_1_420<-p_rt("B.1.420","COVID_BICEPS_B_1_420","magenta",2020.25,2020.58,"B.1.420")
C37 <- p_rt("C.37","COVID_BICEPS_C_37","#3C5488FF",2021.24,2020.55,"Lambda")
B_1_1_7 <- p_rt("B.1.1.7","COVID_BICEPS_B_1_1_7","#DC0000FF",2021.12,2020.67,"Alpha")
P1 <-p_rt("P.1","COVID_BICEPS_P_1","#00A087FF",2021.01,2020.75,"Gamma")
B_1_621 <-p_rt("B.1.621","COVID_BICEPS_B_1_621","#F39B7FFF",2021.03,2020.78,"Mu")
Delta<- p_rt("Delta","COVID_BICEPS_Delta","#4DBBD5B2",2020.93,2021.16,"Delta")
Omicron <- p_rt("Omicron","COVID_BICEPS_Omicron","#8491B4FF",2021.89,2021.69,"Omicron")



# Variables that could affect SARS-CoV-2 transmission 


setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Database/metadata")

colombia2020 <- as.data.frame(fread("data_download_file_reference_2020.csv"))
colombia2021 <- as.data.frame(fread("data_download_file_reference_2021.csv"))
colombia2022 <- as.data.frame(fread("data_download_file_reference_2022.csv"))

colombia <- rbind(colombia2020,colombia2021,colombia2022)

colombia <- colombia %>% filter(location_name == "Colombia")


table(colombia$cumulative_all_fully_vaccinated)

mobility <- ggplot(colombia, aes(date, mobility_mean)) +
  geom_line() +
  theme_classic()+
  #scale_fill_brewer(palette="Spectral")+
  scale_color_npg()+
  scale_fill_npg()+
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 18), axis.title.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", size=18))+
  geom_vline(xintercept =2020-08-01, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  geom_vline(xintercept =2021-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  geom_vline(xintercept =2021-06-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  geom_vline(xintercept =2022-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  xlab("Date") +
  ylab ("Mobility") 


mask <- ggplot(colombia, aes(date, mask_use_mean)) +
  geom_line() +
  theme_classic()+
  #scale_fill_brewer(palette="Spectral")+
  scale_color_npg()+
  scale_fill_npg()+
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 18), axis.title.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", size=18))+
  geom_vline(xintercept =2020-08-01, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  geom_vline(xintercept =2021-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  geom_vline(xintercept =2021-06-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  geom_vline(xintercept =2022-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  xlab("Date") +
  ylab ("Mask use") 


vaccinated <- ggplot(colombia, aes(date, cumulative_all_vaccinated)) +
  geom_step(color="green")+
  geom_line() +
  geom_line(aes(x=date,y=cumulative_all_fully_vaccinated), lty=2)+
  theme_classic()+
  #scale_fill_brewer(palette="Spectral")+
  scale_color_npg()+
  scale_fill_npg()+
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 18), axis.title.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", size=18))+
  geom_vline(xintercept =2020-08-01, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  geom_vline(xintercept =2021-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  geom_vline(xintercept =2021-06-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  geom_vline(xintercept =2022-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  xlab("Date") +
  ylab ("% Vaccinated") 

#####lockown index 


lockdown <- as.data.frame(fread("covid-stringency-index.csv"))

lockdown <- lockdown %>% filter(Entity == "Colombia")

stringency_index <- ggplot(lockdown, aes(Day, stringency_index)) +
  geom_line(color="green", lty=1) +
  theme_classic()+
  #scale_fill_brewer(palette="Spectral")+
  scale_color_npg()+
  scale_fill_npg()+
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 15), axis.title.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", size=18))+
  geom_vline(xintercept =2020-08-01, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  geom_vline(xintercept =2021-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  geom_vline(xintercept =2021-06-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  geom_vline(xintercept =2022-01-15, colour = "grey", lty=1, size=3, alpha=0.5)+
  xlab("Date") +
  ylab ("COVID-19 Stringency Index") 



setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/Figure5")
pdf ('Figure5.pdf',width=25,height=10)
#par(mfrow=c(2,4))

ggarrange(B_1,B_1_111,B_1_1_348,B_1_420,C37,B_1_1_7, B_1_621, P1,Delta,Omicron,mobility,mask,vaccinated,stringency_index, ncol = 4, nrow = 4)

dev.off()


ggsave("Figure4.jpg",plot=p_all,dpi=500,width=30,height=30)


