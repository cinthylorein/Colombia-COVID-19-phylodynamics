library(ggtree)
library(ggplot2)
library(ggimage)
library("igraph")
library(coda)
library(scales)
library(reshape2)
library("scales")
library(ggpubr)
library(extrafont)
library(lubridate)
library(dplyr)

#lineages earliest submission date according to GISAID
B1ed <- 2020.3675
B1111ed <- 2020.3811
B1420ed <- 2020.3784
B11348ed <- 2021.0205
alphaed <- 2021.2918
gammaed <- 2021.0808
mued <- 2021.1904
lambdaed <- 2021.3795
deltaed <- 2021.5603
omicroned <- 2021.9712



#BICEPS

Biceps_B1 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/B1/BICEPS.log", header=TRUE, sep="\t")
Biceps_B_1_111 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/B111/BICEPS.log", header=TRUE, sep="\t")
Biceps_B_1_1_348 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/B11348/BICEPS.log", header=TRUE, sep="\t")
Biceps_B_1_420 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/B1420/BICEPS.log", header=TRUE, sep="\t")
Biceps_B_1_1_7 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/Alpha/BICEPS.log", header=TRUE, sep="\t")
Biceps_B_1_621 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/Mu/BICEPS.log", header=TRUE, sep="\t")
Biceps_C_37 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/Lambda/BICEPS.log", header=TRUE, sep="\t")
Biceps_P_1 <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/Gamma/BICEPS.log", header=TRUE, sep="\t")
Biceps_Delta <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/Delta/BICEPS.log", header=TRUE, sep="\t")
Biceps_Omicron <- read.table(file="~/Documents/Nature Manuscript/BICEPS Logs/Omicron/BICEPS.log", header=TRUE, sep="\t")

length(Biceps_B1)

Biceps_B1 = Biceps_B1[1:2080,]
Biceps_B_1_111 = Biceps_B_1_111[1:2080,]
Biceps_B_1_420 = Biceps_B_1_420[1:2080,]
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
B_1_420 <- as.data.frame(Biceps_B_1_420$TreeHeight)
B_1_420$ID <- factor(c("B.1.420"))
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








####### TMCRA for each variant 

A <- ggplot(NULL) +
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
  geom_vline(xintercept =2020.021, colour = "#91D1C2FF", lty=2)+ #B.1
  geom_vline(xintercept =2020.2, colour = "#64CC80", lty=2)+ #B.1.111
  geom_vline(xintercept =2020.33, colour = "#80A680", lty=2)+ #B.1.1.348
  geom_vline(xintercept =2021.12, colour = "#DC0000FF", lty=2)+#Alpha 
  geom_vline(xintercept =2020.95, colour = "#F39B7FFF", lty=2)+#Mu
  geom_vline(xintercept =2021.24, colour = "#3C5488FF", lty=2)+ #Lambda
  geom_vline(xintercept =2021.01, colour = "#00A087FF", lty=2)+#Gamma
  geom_vline(xintercept =2020.93, colour = "#4DBBD5B2", lty=2)+#Delta
  geom_vline(xintercept =2021.89, colour = "#8491B4FF", lty=2)+#Omicron
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#64CC80",
                                    "B.1.420" = "magenta",
                                    "B.1.1.348" = "#80A680",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
  scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                  "B.1.111" = "#64CC80",
                                  "B.1.420" = "magenta",
                                  "B.1.1.348" = "#80A680",
                                  "Alpha" = "#DC0000FF",
                                  "Mu" = "#F39B7FFF",
                                  "Lambda" = "#3C5488FF",
                                  "Gamma" = "#00A087FF",
                                  "Delta" = "#4DBBD5B2",
                                  "Omicron" = "#8491B4FF"))+
  theme_set(theme_classic(base_size=12))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), plot.title = element_text(face="bold", size=18),
          text = element_text(family = "Times New Roman"))
A

###### Lag : difference between frist report and TMRCA

library(ggridges)

B <- ggplot(NULL) +
  geom_density_ridges(data=B1, aes(x=(B1ed-Biceps_B1$TreeHeight)*365, y=ID, color="B.1", fill="B.1", scale = 8))+
  geom_density_ridges(data=B_1_111, aes(x=(B1111ed-Biceps_B_1_111$TreeHeight)*365, y=ID, color="B.1.111", fill="B.1.111", scale = 8))+
  geom_density_ridges(data=B_1_1_348, aes(x=(B11348ed-Biceps_B_1_1_348$TreeHeight)*365, y=ID, color="B.1.1.348", fill="B.1.1.348", scale = 8))+
  geom_density_ridges(data=B_1_1_7, aes(x=(alphaed-Biceps_B_1_1_7$TreeHeight)*365, y=ID, color="Alpha", fill="Alpha", scale = 8))+
  geom_density_ridges(data=B_1_621, aes(x=(mued-Biceps_B_1_621$TreeHeight)*365, y=ID, color="Mu", fill="Mu", scale = 8))+
  geom_density_ridges(data=C_37, aes(x=(lambdaed-Biceps_C_37$TreeHeight)*365, y=ID, color="Lambda", fill="Lambda", scale = 8))+
  geom_density_ridges(data=P_1, aes(x=(gammaed-Biceps_P_1$TreeHeight)*365, y=ID, color="Gamma", fill="Gamma", scale = 8))+
  geom_density_ridges(data=Delta, aes(x=(deltaed-Biceps_Delta$TreeHeight)*365, y=ID, color="Delta", fill="Delta", scale=8))+
  geom_density_ridges(data=Omi, aes(x=(omicroned-Biceps_Omicron$TreeHeight)*365, y=ID, color="Omicron", fill="Omicron", scale = 8))+
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
  # geom_vline(xintercept =2020.021, colour = "#91D1C2FF", lty=2)+
  # geom_vline(xintercept =2020.2, colour = "#91D1C2FF", lty=2)+
  # geom_vline(xintercept =2020.33, colour = "#91D1C2FF", lty=2)+
  # geom_vline(xintercept =2021.12, colour = "#DC0000FF", lty=2)+#Alpha
  # geom_vline(xintercept =2020.95, colour = "#F39B7FFF", lty=2)+#Mu
  # geom_vline(xintercept =2021.24, colour = "#3C5488FF", lty=2)+ #Lambda
  # geom_vline(xintercept =2021.01, colour = "#00A087FF", lty=2)+#Gamma
  # geom_vline(xintercept =2020.93, colour = "#4DBBD5B2", lty=2)+#Delta
  # geom_vline(xintercept =2021.89, colour = "#8491B4FF", lty=2)+#Omicron
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#64CC80",
                                    "B.1.1.348" = "#80A680",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
 scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                 "B.1.111" = "#64CC80",
                                 "B.1.1.348" = "#80A680",
                                 "Alpha" = "#DC0000FF",
                                 "Mu" = "#F39B7FFF",
                                 "Lambda" = "#3C5488FF",
                                 "Gamma" = "#00A087FF",
                                 "Delta" = "#4DBBD5B2",
                                 "Omicron" = "#8491B4FF"))+
  theme_set(theme_classic(base_size=12))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), plot.title = element_text(face="bold", size=18),
        text = element_text(family = "Times New Roman"))

B

 
#Plot S4C: time between sampling and submission

GISAID = read.csv("~/Documents/R scripts/Phylodynamics/Multinomial fit databases/GISAID.csv", sep = ",")
GISAID = read.csv("~/Documents/Nature Manuscript/GISAID_DATES.csv", sep = ";")
Sub_date = read.csv("~/Documents/Nature Manuscript/Sub_dates.csv", sep = ";")

dates <- merge(GISAID, Sub_date, by="Accession.ID")
colnames(dates) = c("ID", "EPI", "Collection date", "Location", "Host", "Additional Location", "Sampling Strategy", "Gender", "Age",
                    "Status", "Last Vaccinated", "Passage", "Specimen", "Additional Host Info", "Lineage", "Clade", "AA subs", "Collection Date Y",
                    "Submission date", "location.y")

dates$`Collection date` <- parse_date_time(dates$`Collection date`,c('dmY','dmy'))
dates$`Submission date` <- parse_date_time(dates$`Submission date`,c('dmY','dmy'))

#Let's unify sublineages that belong to a single variant

dates$Lineage[grepl("B.1.621|B.1.621.1|B.1.621.2|BB.2", dates$Lineage)] = "Mu"
dates$Lineage[grepl("P.1|P.1.1|P.1.2",dates$Lineage)] ="Gamma"
dates$Lineage[grepl("B.1.1.7",dates$Lineage)] = "Alpha"
dates$Lineage[grepl("B.1.625", dates$Lineage)]
dates$Lineage[grepl("B.1.617.2|AY.", dates$Lineage)] = "Delta"
dates$Lineage[grepl("C.37",dates$Lineage)] = "Lambda"
dates$Lineage[grepl("B.1.1.529|BA.1|BA.2|BA.1.|BA.2.",dates$Lineage)] = "Omicron"

#Filter lineages of interests

lineages_of_i <- c("B.1", "B.1.111", "B.1.420", "B.1.1.348", "Alpha", "Gamma", "Lambda", "Mu", "Delta", "Omicron")

dates = dates[dates$Lineage %in% lineages_of_i,]

#convert all dates into decimal years and subtract Collection date to submission date

dates$`Collection date` = decimal_date(dates$`Collection date`)
dates$`Submission date` = decimal_date(dates$`Submission date`)

dates$Difference = dates$`Submission date`-dates$`Collection date`
dates$Difference_days = dates$Difference*365

dates %>% group_by(Lineage)

dates$Lineage = factor(dates$Lineage, levels = c("B.1", "B.1.111", "B.1.420", "B.1.1.348", "Alpha", "Gamma", "Mu", "Lambda", "Delta", "Omicron"))

C <- ggplot(dates, aes(x = Difference_days, y = Lineage)) +
  geom_boxplot(width=0.3,aes(fill=Lineage, alpha=0.7),outlier.shape = NA, position =position_dodge(width=1)) +
  geom_point(size = 0.5, alpha = .3, position = position_jitter(seed = 1, width = 0.1, height = 0.1), aes(color = Lineage, fill = Lineage))+ 
  ggdist::stat_halfeye(
    adjust = .5,
    width = 0.4,
    height = 0.5,
    justification = -0.45,
    .width = 0,
    point_colour = NA,
    aes(fill = Lineage, alpha = 0.7))+
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#64CC80",
                                    "B.1.1.348" = "#80A680",
                                    "B.1.420" = "magenta",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
  scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                  "B.1.111" = "#64CC80",
                                  "B.1.1.348" = "#80A680",
                                  "B.1.420" = "magenta",
                                  "Alpha" = "#DC0000FF",
                                  "Mu" = "#F39B7FFF",
                                  "Lambda" = "#3C5488FF",
                                  "Gamma" = "#00A087FF",
                                  "Delta" = "#4DBBD5B2",
                                  "Omicron" = "#8491B4FF"))+
  scale_x_continuous(name = "Mean days to submission(days)") +
  theme_set(theme_classic(base_size=12))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), plot.title = element_text(face="bold", size=18),
        text = element_text(family = "Times New Roman"),
        legend.position = "none")
C

setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/Figures")

p_all<-ggarrange(A,labels = c("", "A"), ncol = 2, ggarrange(B,C, nrow = 2, 
                 labels = c("B", "C"),
                 legend = "none"),
                 font.label = list(size = 14, color = "black", face = "bold", family = "Times New Roman"),
                 common.legend = TRUE,
                 legend = "left")
p_all
ggsave("S4Fig.jpg",plot=p_all,dpi=900,width=25,height=15, units = "cm")
ggsave("S4Fig.svg",plot=p_all,dpi=900,width=15,height=5)
ggsave("S4Fig.pdf",plot=p_all,dpi=900,width=25,height=15, units = "cm")

