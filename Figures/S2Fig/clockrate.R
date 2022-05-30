library(ggtree)
library(ggplot2)
library(ggimage)
library("igraph")
library(coda)
library(scales)
library(reshape2)
library("scales")
library(ggpubr)


#BICEPS

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/B_1")
Biceps_B1 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/B_1_111")
Biceps_B_1_111 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/B_1_1_348")
Biceps_B_1_1_348 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/B_1_420")
Biceps_B_1_420 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/B_1_1_7")
Biceps_B_1_1_7 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/B_1_621")
Biceps_B_1_621 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/C_37")
Biceps_C_37 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/P_1")
Biceps_P_1 <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/Delta")
Biceps_Delta <- read.table(file="BICEPS.log", header=TRUE, sep="\t")
setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Biceps/Omicron")
Biceps_Omicron <- read.table(file="BICEPS.log", header=TRUE, sep="\t")

length(Biceps_B1)

Biceps_B1 = Biceps_B1[1:2080,]
Biceps_B_1_111 = Biceps_B_1_111[1:2080,]
Biceps_B_1_1_348 = Biceps_B_1_1_348[1:2080,]
Biceps_B_1_1_7 = Biceps_B_1_1_7[1:2080,]
Biceps_B_1_621 = Biceps_B_1_621[1:2080,]
Biceps_C_37 = Biceps_C_37[1:2080,]
Biceps_P_1 = Biceps_P_1[1:2080,]
Biceps_Delta = Biceps_Delta[1:2080,]
Biceps_Omicron = Biceps_Omicron[1:2080,]


mrs <- 215/365 + 2021
Biceps_B1$TreeHeight = mrs  - Biceps_B1$TreeHeight
mrs <- 163/365 + 2021
Biceps_B_1_111$TreeHeight = mrs  - Biceps_B_1_111$TreeHeight
mrs <- 138/365 + 2021
Biceps_B_1_1_348$TreeHeight = mrs  - Biceps_B_1_1_348$TreeHeight
mrs <- 213/365 + 2021
Biceps_B_1_1_7$TreeHeight = mrs  - Biceps_B_1_1_7$TreeHeight
mrs <- 329/365 + 2021
Biceps_B_1_621$TreeHeight = mrs  - Biceps_B_1_621$TreeHeight
mrs <- 156/365 + 2021
Biceps_C_37$TreeHeight = mrs  - Biceps_C_37$TreeHeight
mrs <- 332/365 + 2021
Biceps_P_1$TreeHeight = mrs  - Biceps_P_1$TreeHeight
mrs <- 17/365 + 2022
Biceps_Delta$TreeHeight = mrs  - Biceps_Delta$TreeHeight
mrs <- 17/365 + 2022
Biceps_Omicron$TreeHeight = mrs  - Biceps_Omicron$TreeHeight


B1 <- as.data.frame(Biceps_B1$TreeHeight)
B1$ID <- factor(c('B.1'))
B_1_111 <- as.data.frame(Biceps_B_1_111$TreeHeight)
B_1_111$ID <- factor(c('B.1.111'))
B_1_1_348 <- as.data.frame(Biceps_B_1_1_348$TreeHeight)
B_1_1_348$ID <- factor(c('B.1.348'))
B_1_1_7 <- as.data.frame(Biceps_B_1_1_7$TreeHeight)
B_1_1_7$ID <- factor(c('Alpha'))
B_1_621 <- as.data.frame(Biceps_B_1_621$TreeHeight)
B_1_621$ID <- factor(c('Mu'))
C_37 <- as.data.frame(Biceps_C_37$TreeHeight)
C_37$ID <- factor(c('Lambda'))
P_1 <- as.data.frame(Biceps_P_1$TreeHeight)
P_1$ID <- factor(c('Gamma'))
Delta <- as.data.frame(Biceps_Delta$TreeHeight)
Delta$ID <- factor(c('Delta'))
Omi<- as.data.frame(Biceps_Omicron$TreeHeight)
Omi$ID <- factor(c('Omicron'))

Omi_Re <- cbind(Omi_Re, Omi_times)
## Function getting nean and HPD

my_summary = function(s) {
  q = HPDinterval(mcmc(s))
  return (data.frame(y = median(s), ymin = q[1], ymax = q[2]))
}


meanexp_Biceps_B1_TreeHeight<- my_summary(Biceps_B1$TreeHeight)
fr <- 510/365
(meanexp_Biceps_B1_TreeHeight-fr)*365

meanexp_Biceps_B_1_111_TreeHeight <- my_summary(Biceps_B_1_111$TreeHeight)
fr <- 454/365
(meanexp_Biceps_B_1_111_TreeHeight-fr)*365

meanexp_Biceps_B_1_1_348_TreeHeight <- my_summary(Biceps_B_1_1_348$TreeHeight)
fr <- 384/365
(meanexp_Biceps_B_1_1_348_TreeHeight-fr)*365

meanexp_Biceps_B_1_1_7_TreeHeight <- my_summary(Biceps_B_1_1_7$TreeHeight)
fr <- 213/365
(meanexp_Biceps_B_1_1_7_TreeHeight-fr)*365

meanexp_Biceps_B_1_621_TreeHeight <- my_summary(Biceps_B_1_621$TreeHeight)
fr <- 330/365
(meanexp_Biceps_B_1_621_TreeHeight-fr)*365

meanexp_Biceps_C_37_TreeHeight <- my_summary(Biceps_C_37$TreeHeight)
fr <- 156/365
(meanexp_Biceps_C_37_TreeHeight-fr)*365

meanexp_Biceps_P_1_TreeHeight <- my_summary(Biceps_P_1$TreeHeight)
fr <- 332/365
(meanexp_Biceps_P_1_TreeHeight-fr)*365

meanexp_Biceps_Delta_TreeHeight <- my_summary(Biceps_Delta$TreeHeight)
fr <- 252/365
(meanexp_Biceps_Delta_TreeHeight-fr)*365

meanexp_Biceps_Omicron_TreeHeight <- my_summary(Biceps_Omicron$TreeHeight)
fr <- 42/365
(meanexp_Biceps_Omicron_TreeHeight-fr)*365


setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Figures/clock_rate")

p_all<-ggarrange(sky,BDSK,BICEPS, ncol = 3, nrow = 1)
ggsave("S2Fig.jpg",plot=p_all,dpi=500,width=15,height=5)





####### TMCRA for each variant 

BICEPS <- ggplot(NULL) +
  geom_density(data=Biceps_B1, aes(TreeHeight, color="B.1" ))+
  geom_density(data=Biceps_B_1_111, aes(TreeHeight, color="B.1.111"))+
  geom_density(data=Biceps_B_1_1_348, aes(TreeHeight, color="B.1.1.348"))+
  geom_density(data=Biceps_B_1_1_7, aes(TreeHeight, color="Alpha"))+
  geom_density(data=Biceps_B_1_621, aes(TreeHeight, color="Mu"))+
  geom_density(data=Biceps_C_37, aes(TreeHeight, color="Lambda"))+
  geom_density(data=Biceps_P_1, aes(TreeHeight, color="Gamma"))+
  geom_density(data=Biceps_Delta, aes(TreeHeight, color="Delta"))+
  geom_density(data=Biceps_Omicron, aes(TreeHeight, color="Omicron"))+
  scale_x_continuous(name = "TMRCA") +
  #scale_y_continuous(name = "Density")+
  #labs(title = "BICEPS model")+
  #wave1 Aug 2020 214/365
  geom_vline(xintercept =2020.58, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  geom_vline(xintercept =2021.06, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  geom_vline(xintercept =2021.41, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  geom_vline(xintercept =2022.005, colour = "grey", lty=1, size=3, alpha=0.5)+
  #First colombian case reported 
  geom_vline(xintercept =2020.021, colour = "#91D1C2FF", lty=2)+
  geom_vline(xintercept =2020.2, colour = "#91D1C2FF", lty=2)+
  geom_vline(xintercept =2020.33, colour = "#91D1C2FF", lty=2)+ 
  geom_vline(xintercept =2021.12, colour = "#DC0000FF", lty=2)+#Alpha 
  geom_vline(xintercept =2020.95, colour = "#F39B7FFF", lty=2)+#Mu
  geom_vline(xintercept =2021.24, colour = "#3C5488FF", lty=2)+ #Lambda
  geom_vline(xintercept =2021.01, colour = "#00A087FF", lty=2)+#Gamma
  geom_vline(xintercept =2020.93, colour = "#4DBBD5B2", lty=2)+#Delta
  geom_vline(xintercept =2021.89, colour = "#8491B4FF", lty=2)+#Omicron
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#91D1C2FF",
                                    "B.1.1.348" = "#91D1C2FF",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
  scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                  "B.1.111" = "#91D1C2FF",
                                  "B.1.1.348" = "#91D1C2FF",
                                  "Alpha" = "#DC0000FF",
                                  "Mu" = "#F39B7FFF",
                                  "Lambda" = "#3C5488FF",
                                  "Gamma" = "#00A087FF",
                                  "Delta" = "#4DBBD5B2",
                                  "Omicron" = "#8491B4FF"))+
  theme_set(theme_classic(base_size=12))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), plot.title = element_text(face="bold", size=18))

###### Lag : difference between frist report and TMRCA

library(ggridges)

BICEPS <- ggplot(NULL) +
  geom_density_ridges(data=B1, aes(x=Biceps_B1$TreeHeight*365-510, y=ID, color="B.1", fill="B.1", scale = 8))+
  geom_density_ridges(data=B_1_111, aes(x=Biceps_B_1_111$TreeHeight*365-454, y=ID, color="B.1.111", fill="B.1.111", scale = 8))+
  geom_density_ridges(data=B_1_1_348, aes(x=Biceps_B_1_1_348$TreeHeight*365-384, y=ID, color="B.1.1.348", fill="B.1.1.348", scale = 8))+
  geom_density_ridges(data=B_1_1_7, aes(x=Biceps_B_1_1_7$TreeHeight*365-(213+30), y=ID, color="Alpha", fill="Alpha", scale = 8))+
  geom_density_ridges(data=B_1_621, aes(x=Biceps_B_1_621$TreeHeight*365-330, y=ID, color="Mu", fill="Mu", scale = 8))+
  geom_density_ridges(data=C_37, aes(x=Biceps_C_37$TreeHeight*365-(156+20), y=ID, color="Lambda", fill="Lambda", scale = 8))+
  geom_density_ridges(data=P_1, aes(x=Biceps_P_1$TreeHeight*365-(332+10), y=ID, color="Gamma", fill="Gamma", scale = 8))+
  #geom_density_ridges(data=Delta, aes(y=Biceps_Delta$TreeHeight*365-(252+60), x=ID, color="Delta"))+
  geom_density_ridges(data=Omi, aes(x=Biceps_Omicron$TreeHeight*365-(42+10), y=ID, color="Omicron", fill="Omicron", scale = 8))+
  scale_x_continuous(name = "Detection lag(days)") +
  #scale_y_continuous(name = "Detection lag(days)")+
  #labs(title = "BICEPS model")+
  #wave1 Aug 2020 214/365
  #geom_vline(xintercept =2020.58, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  #geom_vline(xintercept =2021.06, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  #geom_vline(xintercept =2021.41, colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  #geom_vline(xintercept =2022.005, colour = "grey", lty=1, size=3, alpha=0.5)+
  #First colombian case reported 
  #geom_vline(xintercept =2020.021, colour = "#91D1C2FF", lty=2)+
  #geom_vline(xintercept =2020.2, colour = "#91D1C2FF", lty=2)+
  #geom_vline(xintercept =2020.33, colour = "#91D1C2FF", lty=2)+ 
  #geom_vline(xintercept =2021.12, colour = "#DC0000FF", lty=2)+#Alpha 
  #geom_vline(xintercept =2020.95, colour = "#F39B7FFF", lty=2)+#Mu
  #geom_vline(xintercept =2021.24, colour = "#3C5488FF", lty=2)+ #Lambda
  #geom_vline(xintercept =2021.01, colour = "#00A087FF", lty=2)+#Gamma
  #geom_vline(xintercept =2020.93, colour = "#4DBBD5B2", lty=2)+#Delta
  #geom_vline(xintercept =2021.89, colour = "#8491B4FF", lty=2)+#Omicron
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#91D1C2FF",
                                    "B.1.1.348" = "#91D1C2FF",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
 scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                  "B.1.111" = "#91D1C2FF",
                                 "B.1.1.348" = "#91D1C2FF",
                                  "Alpha" = "#DC0000FF",
                                 "Mu" = "#F39B7FFF",
                                  "Lambda" = "#3C5488FF",
                                  "Gamma" = "#00A087FF",
                                  "Delta" = "#4DBBD5B2",
                                 "Omicron" = "#8491B4FF"))+
  theme_set(theme_classic(base_size=12))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), plot.title = element_text(face="bold", size=18))

