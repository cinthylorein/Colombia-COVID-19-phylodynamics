# BDSKY analysis
# https://taming-the-beast.org/tutorials/Skyline-plots/

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#BiocManager::install(bdskytools)


if (!require(devtools, quietly = T)) install.packages("devtools")
install.packages("remotes")
remotes::install_github("laduplessis/bdskytools")
library(devtools)
library(bdskytools)
library(seqinr)
require(data.table)
library(ggplot2)
library(ggpubr)
library("scales")


#### B_1

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1/10_dimensions")
fname <- "BDSKL2_COMB.log"    
B_1_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_Re_sky    <- getSkylineSubset(B_1_lf, "reproductiveNumber")
B_1_Re_hpd    <- getMatrixHPD(B_1_Re_sky)
B_1_delta_hpd <- getHPD(B_1_lf$becomeUninfectiousRate)

# plot the smooth skyline 510/730
period = 1.55
timegrid       <- seq(0,period,length.out=72)
B_1_Re_gridded     <- gridSkyline(B_1_Re_sky, B_1_lf$origin, timegrid)
B_1_Re_gridded_hpd <- getMatrixHPD(B_1_Re_gridded)

# last sample date 2021-08-03 (215)/365
date = 215/365 + 2021
B_1_times <- date - timegrid


#### B_1_1_111

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_111/10_dimensions")
fname <- "BDSKL2_COMB.log"    
B_1_1_111_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_1_111_Re_sky    <- getSkylineSubset(B_1_1_111_lf, "reproductiveNumber")
B_1_1_111_Re_hpd    <- getMatrixHPD(B_1_1_111_Re_sky)
B_1_1_111_delta_hpd <- getHPD(B_1_1_111_lf$becomeUninfectiousRate)

# plot the smooth skyline 454/730
period = 1.28
timegrid       <- seq(0,period,length.out=64)
B_1_1_111_Re_gridded     <- gridSkyline(B_1_1_111_Re_sky, B_1_1_111_lf$origin, timegrid)
B_1_1_111_Re_gridded_hpd <- getMatrixHPD(B_1_1_111_Re_gridded)

# last sample date 2021-06-12 (163)/365
date = 163/365 + 2021
B_1_1_111_times <- date - timegrid

#### B_1_1_348

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_1_348/10_dimensions")
fname <- "BDSKL2_COMB.log"    
B_1_1_348_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_1_348_Re_sky    <- getSkylineSubset(B_1_1_348_lf, "reproductiveNumber")
B_1_1_348_Re_hpd    <- getMatrixHPD(B_1_1_348_Re_sky)
B_1_1_348_delta_hpd <- getHPD(B_1_1_348_lf$becomeUninfectiousRate)

# plot the smooth skyline 384/730
period = 1.11
timegrid       <- seq(0,period,length.out=54)
B_1_1_348_Re_gridded     <- gridSkyline(B_1_1_348_Re_sky, B_1_1_348_lf$origin, timegrid)
B_1_1_348_Re_gridded_hpd <- getMatrixHPD(B_1_1_348_Re_gridded)

# last sample date 2021-05-18 (138)/365
date = 138/365 + 2021
B_1_1_348_times <- date - timegrid

#### B_1_420

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_420")
fname <- "BDSKL2_COMB.log"    
B_1_420_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_420_Re_sky    <- getSkylineSubset(B_1_420_lf, "reproductiveNumber")
B_1_420_Re_hpd    <- getMatrixHPD(B_1_420_Re_sky)
B_1_420_delta_hpd <- getHPD(B_1_420_lf$becomeUninfectiousRate)

# plot the smooth skyline 118/365
period = 0.32
timegrid       <- seq(0,period,length.out=16)
B_1_420_Re_gridded     <- gridSkyline(B_1_420_Re_sky, B_1_420_lf$origin, timegrid)
B_1_420_Re_gridded_hpd <- getMatrixHPD(B_1_420_Re_gridded)

# last sample date 2021-05-18 (138)/365
date = 118/365 + 2021
B_1_420_times <- date - timegrid


#### Delta

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/Delta/10_dimensions")
fname <- "BDSKL2_COMB.log"    
Delta_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
Delta_Re_sky    <- getSkylineSubset(Delta_lf, "reproductiveNumber")
Delta_Re_hpd    <- getMatrixHPD(Delta_Re_sky)
Delta_delta_hpd <- getHPD(Delta_lf$becomeUninfectiousRate)

# plot the smooth skyline 407/1095
period = 1.68
timegrid       <- seq(0,period,length.out=58)
Delta_Re_gridded     <- gridSkyline(Delta_Re_sky, Delta_lf$origin, timegrid)
Delta_Re_gridded_hpd <- getMatrixHPD(Delta_Re_gridded)

# last sample date 2021-06-10 (189)/365
date = 17/365 + 2022
Delta_times <- date - timegrid


#P.1

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/P_1/10_dimensions")
fname <- "BDSKL2_COMB.log"    
lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
P1_Re_sky    <- getSkylineSubset(lf, "reproductiveNumber")
P1_Re_hpd    <- getMatrixHPD(P1_Re_sky)
P1_hpd <- getHPD(lf$becomeUninfectiousRate)


# plot the smooth skyline 332/365
period = 1.07
timegrid       <- seq(0,period,length.out=47)
P1_Re_gridded     <- gridSkyline(P1_Re_sky,    lf$origin, timegrid)
P1_Re_gridded_hpd <- getMatrixHPD(P1_Re_gridded)

# last sample date 
date = 332/365 + 2021
P1_times <- date - timegrid

#### Alpha

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_1_7/10_dimensions")
fname <- "BDSKL2_COMB.log"    
B_1_1_7_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_1_7_Re_sky    <- getSkylineSubset(B_1_1_7_lf, "reproductiveNumber")
B_1_1_7_Re_hpd    <- getMatrixHPD(B_1_1_7_Re_sky)
B_1_1_7_delta_hpd <- getHPD(B_1_1_7_lf$becomeUninfectiousRate)

# plot the smooth skyline 213/365
period = 0.87
timegrid       <- seq(0,period,length.out=30)
B_1_1_7_Re_gridded     <- gridSkyline(B_1_1_7_Re_sky, B_1_1_7_lf$origin, timegrid)
B_1_1_7_Re_gridded_hpd <- getMatrixHPD(B_1_1_7_Re_gridded)

# last sample date 
date = 213/365 + 2021
B_1_1_7_times <- date - timegrid



#### C37

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/C_37/10_dimensions")
fname <- "BDSKL2_COMB.log"    
C37_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
C37_Re_sky    <- getSkylineSubset(C37_lf, "reproductiveNumber")
C37_Re_hpd    <- getMatrixHPD(C37_Re_sky)
C37_delta_hpd <- getHPD(C37_lf$becomeUninfectiousRate)

# plot the smooth skyline 156/365
period = 0.73
timegrid       <- seq(0,period,length.out=22)
C37_Re_gridded     <- gridSkyline(C37_Re_sky,    C37_lf$origin, timegrid)
C37_Re_gridded_hpd <- getMatrixHPD(C37_Re_gridded)

# last sample date
date = 156/365 + 2021
C37_times <- date - timegrid


#### Omicron

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/Omicron/10_dimensions")
fname <- "BDSKL2_COMB.log"    
Omi_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
Omi_Re_sky    <- getSkylineSubset(Omi_lf, "reproductiveNumber")
Omi_Re_hpd    <- getMatrixHPD(Omi_Re_sky)
Omi_delta_hpd <- getHPD(Omi_lf$becomeUninfectiousRate)

plotSkyline(1:10, Omi_Re_hpd, type='step', ylab="R")

# plot the smooth skyline 42/730
period = 0.2
timegrid       <- seq(0,period,length.out=12)
Omi_Re_gridded     <- gridSkyline(Omi_Re_sky, Omi_lf$origin, timegrid)
Omi_Re_gridded_hpd <- getMatrixHPD(Omi_Re_gridded)

# last sample date 2021-06-10 (189)/365
date = 17/365 + 2022
Omi_times <- date - timegrid


#### SP

setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/SP")
fname <- "BDSKL2_COMB.log"    
SP_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
SP_Re_sky    <- getSkylineSubset(SP_lf, "reproductiveNumber")
SP_Re_hpd    <- getMatrixHPD(SP_Re_sky)
SP_delta_hpd <- getHPD(SP_lf$becomeUninfectiousRate)

# plot the smooth skyline
period = 3.83
timegrid       <- seq(0,period,length.out=34)
SP_Re_gridded     <- gridSkyline(SP_Re_sky,    SP_lf$origin, timegrid)
SP_Re_gridded_hpd <- getMatrixHPD(SP_Re_gridded)

# last sample date 2021-06-10 (189)/365
date = 233/365 + 2021
SP_times <- date - timegrid





setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures/figure4")

############################# Plot Figure 3 raw points  


############## Re plot smothing 

setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Figures")
pdf ('Figure4.pdf',width=20,height=10)
par(mfrow=c(3,4))


####1.B_1



plotSkylinePretty(B_1_times, B_1_Re_gridded_hpd, type='smooth',  main="B.1", axispadding=0.0, 
                  col=pal.dark(cgreen, 0.6), fill=pal.dark(cgreen, 0.1), 
                  ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
#line(data= Mob, aes(x=time,y=mobility_mean), lty=2)
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first squence reported 2020-03-12 (31+29+12)/365=0.19
abline(v=2020.19, col ="red", lty=1, main="B1")
# Lineage first squence reported 2020-01-08 08/365=0.021
abline(v=2020.021, col ="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)


#### 3. B_1_1_111

plotSkylinePretty(B_1_1_111_times, B_1_1_111_Re_gridded_hpd, type='smooth', main="B.1.111", axispadding=0.0, 
                  col=pal.dark(cgreen, 0.6), fill=pal.dark(cgreen, 0.1), 
                  side=2, yline=2.5, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2020-03-13 (31+29+13)/365=0.2
abline(v=2020.2, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2020-03-07 (31+29+07)/365=0.18
abline(v=2020.18, col="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)


#### 5. B_1_1_348


plotSkylinePretty(B_1_1_348_times, B_1_1_348_Re_gridded_hpd, type='smooth', main="B.1.1.348", axispadding=0.0, 
                  col=pal.dark(cgreen, 0.6), fill=pal.dark(cgreen, 0.1), 
                  side=2, yline=2.5, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2020-04-30 (31+29+31+30)/365=0.33
abline(v=2020.33, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2020-04-30 (31+29+31+30)/365=0.33
abline(v=2020.33, col="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)

#### 5. B_1_420


plotSkylinePretty(B_1_420_times, B_1_420_Re_gridded_hpd, type='smooth', main="B.1.420", axispadding=0.0, 
                  col=pal.dark(cred, 0.6), fill=pal.dark(cred, 0.1), 
                  side=2, yline=2.5, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2020-03-16 (31+29+31+30)/365=0.33
abline(v=2020.25, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2020-07-13 (31+29+31+30)/365=0.33
abline(v=2020.58, col="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)


####8.C37

plotSkylinePretty(C37_times, C37_Re_gridded_hpd, type='smooth', main="Lambda", axispadding=0.0, 
                  col=pal.dark(cblue, 0.9), fill=pal.dark(cblue, 0.4), 
                  side=2, yline=2.5, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2021-03-30 (31+29+30)/365=0.24
abline(v=2021.24, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2020-07-21 (31+29+31+30+31+30+21)/365=0.55
abline(v=2020.55, col="red", lty=2, main="B1")
#text(mean(x1),(b*max(x1)/2+a),"text",srt=0.2,pos=3)
text(mean(B_1_times),"Lockdown",srt=0.2)


#### 6.B_1_1_7


plotSkylinePretty(B_1_1_7_times, B_1_1_7_Re_gridded_hpd, type='smooth', main="Alpha", axispadding=0.0, 
                  col=pal.dark(cred, 0.9), fill=pal.dark(cred, 0.1), 
                  xlab="Date", ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=FALSE,  
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2021-02-15 46/365=0.12
abline(v=2021.12, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2020-09-03 (31+29+31+30+31+30+31+31+03)/365=0.67
abline(v=2020.67, col="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)

#### 2.P1

plotSkylinePretty(P1_times, P1_Re_gridded_hpd, type='smooth',  main="Gamma", axispadding=0.0, 
                  col=pal.dark(cgreen, 1.0), fill=pal.dark(cgreen, 0.4), 
                  xlab="Date", side=2, yline=2.5, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first squence reported 2021-01-04 (4)/365=0.01
abline(v=2021.01, col ="red", lty=1, main="B1")
# Lineage's first squence reported 2020-10-01 (31+29+31+30+31+30+31+31+30+01)/365=0.75
abline(v=2020.75, col ="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)

#### 7. Delta


plotSkylinePretty(Delta_times, Delta_Re_gridded_hpd, type='smooth', main="Delta", axispadding=0.0,
                  col=pal.dark(cblue, 0.9), fill=pal.dark(cblue, 0.4), 
                  xlab="Time", side=2, yline=2, xline=2, xgrid=FALSE, 
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2020,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2020-12-07 342/365=0.93
abline(v=2020.93, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2021-03-01 60/365=0.16
abline(v=2021.16, col="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)

#### 7. Omicron


plotSkylinePretty(Omi_times, Omi_Re_gridded_hpd, type='smooth', main="Omicron", axispadding=0.0, 
                  col=pal.dark(cpurple, 0.9), fill=pal.dark(cpurple, 0.4), 
                  xlab="Time", side=2, yline=2, xline=2, xgrid=FALSE,  
                  ygrid=FALSE, gridcol=pal.dark(cgray),cex.label = 2.0, cex.axis = 2.0, cex.main=2,
                  ylims=c(0,3), xlims=c(2021.5,2022.1))
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+29+25)/365=0.23
abline(v=2020.23, col="blue", lty=2, main="B1")
# Colombia first case 2020-02-26 (31+26)/365=0.15
abline(v=2020.15, col="black", lty=1, main="B1")
# Colombia first sequence reported 2021-11-22 326/365=0.89
abline(v=2021.89, col="red", lty=1, main="B1")
# Lineage's first sequence reported 2020-09-11 254/365=0.69
abline(v=2021.69, col="red", lty=2, main="B1")
text(mean(B_1_times),"Lockdown",srt=0.2)

dev.off()




setwd("/Volumes/GoogleDrive/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia")
pdf ('Figure_Raw.pdf',width=12,height=15)
par(mfrow=c(5,2))

##1.SP 
plotSkylinePretty(0:10, SP_Re_hpd, type='smooth', ylab="Re", main="Sampling from the prior",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

####2.B_1

plotSkylinePretty(0:10, B_1_Re_hpd, type='smooth', ylab="Re",  main="B1")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

#### 3. B_1_1_348

plotSkylinePretty(0:10, B_1_1_348_Re_hpd, type='smooth', ylab="Re",  main="B_1_1_348")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

#### 4. B_1_1_111

plotSkylinePretty(0:10, B_1_1_111_Re_hpd, type='smooth', ylab="Re",  main="B_1_1_111")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)


#### 5.B_1_625

plotSkylinePretty(0:10, B_1_625_Re_hpd, type='smooth', ylab="Re",  main="B_1_625")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

#### 6.P1

plotSkylinePretty(0:10, P1_Re_hpd, type='smooth', ylab="Re",  main="P1")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

#### 7.B_1_1_7

plotSkylinePretty(0:10, B_1_1_7_Re_hpd, type='smooth', ylab="Re",  main="B_1_1_7")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

#### 8.B_1_621

plotSkylinePretty(0:10, B_1_621_Re_hpd, type='smooth', ylab="Re",  main="B_1_621")
abline(h=1, col="red", lty=2)
#Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

####5.C37

plotSkylinePretty(0:10, C37_Re_hpd, type='smooth', ylab="Re",  main="C37")
abline(h=1, col="red", lty=2)
# Colombia lockdown 2020-03-25 (31+23)/365=0.1479452
#lockdown = 2
#abline(v=lockdown, col="yellow", lty=2)

dev.off()

