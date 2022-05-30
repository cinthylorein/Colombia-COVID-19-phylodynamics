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


library(lubridate)
B_1_times <- format(date_decimal(B_1_times), "%Y-%m-%d")

B_1_Re <- as.data.frame(t(B_1_Re_gridded_hpd))
B_1_Re <- cbind(B_1_Re, B_1_times)
names(B_1_Re) <-c("Lower Bound","mean","Upper Bound","date")


B_1_Re <- as.data.frame(t(B_1_Re_hpd))
B_1_Re$ID <- factor(c('B.1'))

B1Re_10 <- as.data.frame(cbind(B_1_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
              B_1_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.6,
              B_1_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.9,
              B_1_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B1Re_10$ID <- factor(c('B.1'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1/15_dimensions")
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


library(lubridate)
B_1_times <- format(date_decimal(B_1_times), "%Y-%m-%d")

B_1_Re <- as.data.frame(t(B_1_Re_gridded_hpd))
B_1_Re <- cbind(B_1_Re, B_1_times)
names(B_1_Re) <-c("Lower Bound","mean","Upper Bound","date")


B_1_Re <- as.data.frame(t(B_1_Re_hpd))
B_1_Re$ID <- factor(c('B.1'))

B1Re_15 <- as.data.frame(cbind(B_1_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                            B_1_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                            B_1_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                            B_1_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B1Re_15$ID <- factor(c('B.1'))
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

### dataframe with Re values and hpds 

B_1_1_111_times <- format(date_decimal(B_1_1_111_times), "%Y-%m-%d")

B_1_1_111_Re <- as.data.frame(t(B_1_1_111_Re_gridded_hpd))
B_1_1_111_Re <- cbind(B_1_1_111_Re, B_1_1_111_times)
names(B_1_1_111_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_1_111_Re <- as.data.frame(t(B_1_1_111_Re_hpd))
B_1_1_111_Re$ID <- factor(c('B.1.111'))

B111Re_10 <- as.data.frame(cbind(B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B111Re_10$ID <- factor(c('B.1.111'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_111/15_dimensions")
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

### dataframe with Re values and hpds 

B_1_1_111_times <- format(date_decimal(B_1_1_111_times), "%Y-%m-%d")

B_1_1_111_Re <- as.data.frame(t(B_1_1_111_Re_gridded_hpd))
B_1_1_111_Re <- cbind(B_1_1_111_Re, B_1_1_111_times)
names(B_1_1_111_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_1_111_Re <- as.data.frame(t(B_1_1_111_Re_hpd))
B_1_1_111_Re$ID <- factor(c('B.1.111'))

B111Re_15 <- as.data.frame(cbind(B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                 B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                 B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                 B_1_1_111_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B111Re_15$ID <- factor(c('B.1.111'))
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

### dataframe with Re values and hpds 

B_1_1_348_times <- format(date_decimal(B_1_1_348_times), "%Y-%m-%d")

B_1_1_348_Re <- as.data.frame(t(B_1_1_348_Re_gridded_hpd))
B_1_1_348_Re <- cbind(B_1_1_348_Re, B_1_1_348_times)
names(B_1_1_348_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_1_348_Re <- as.data.frame(t(B_1_1_348_Re_hpd))
B_1_1_348_Re$ID <- factor(c('B.1.1.348'))

B11348Re_10 <- as.data.frame(cbind(B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                  B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                  B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                  B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B11348Re_10$ID <- factor(c('B.1.1.348'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_1_348/15_dimensions")
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

### dataframe with Re values and hpds 

B_1_1_348_times <- format(date_decimal(B_1_1_348_times), "%Y-%m-%d")

B_1_1_348_Re <- as.data.frame(t(B_1_1_348_Re_gridded_hpd))
B_1_1_348_Re <- cbind(B_1_1_348_Re, B_1_1_348_times)
names(B_1_1_348_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_1_348_Re <- as.data.frame(t(B_1_1_348_Re_hpd))
B_1_1_348_Re$ID <- factor(c('B.1.1.348'))

B11348Re_15 <- as.data.frame(cbind(B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                   B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                   B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                   B_1_1_348_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B11348Re_15$ID <- factor(c('B.1.1.348'))


#### B_1_420

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_420")
fname <- "BDSKL2_COMB.log"    
B_1_420_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_420_Re_sky    <- getSkylineSubset(B_1_420_lf, "reproductiveNumber")
B_1_420_Re_hpd    <- getMatrixHPD(B_1_420_Re_sky)
B_1_420_delta_hpd <- getHPD(B_1_420_lf$becomeUninfectiousRate)

# plot the smooth skyline 384/730
period = 0.32
timegrid       <- seq(0,period,length.out=16)
B_1_420_Re_gridded     <- gridSkyline(B_1_420_Re_sky, B_1_420_lf$origin, timegrid)
B_1_420_Re_gridded_hpd <- getMatrixHPD(B_1_420_Re_gridded)

# last sample date 2021-05-18 (138)/365
date = 118/365 + 2021
B_1_420_times <- date - timegrid

### dataframe with Re values and hpds 

B_1_420_times <- format(date_decimal(B_1_420_times), "%Y-%m-%d")

B_1_420_Re <- as.data.frame(t(B_1_420_Re_gridded_hpd))
B_1_420_Re <- cbind(B_1_420_Re, B_1_420_times)
names(B_1_420_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_420_Re <- as.data.frame(t(B_1_420_Re_hpd))
B_1_420_Re$ID <- factor(c('B_1_420'))

B_1_420_Re <- as.data.frame(cbind(B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                  B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                  B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                  B_1_420_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B_1_420_Re$ID <- factor(c('B.1.420'))


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

### dataframe with Re values and hpds 

Delta_times <- format(date_decimal(Delta_times), "%Y-%m-%d")

Delta_Re <- as.data.frame(t(Delta_Re_gridded_hpd))
Delta_Re <- cbind(Delta_Re, Delta_times)
names(Delta_Re) <-c("Lower Bound","mean","Upper Bound","date")

Delta_Re <- as.data.frame(t(Delta_Re_hpd))
Delta_Re$ID <- factor(c('Delta'))

DeltaRe_10 <- as.data.frame(cbind(Delta_Re_sky$reproductiveNumber_BDSKY_Serial.1, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.2, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                  Delta_Re_sky$reproductiveNumber_BDSKY_Serial.4, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.5, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                  Delta_Re_sky$reproductiveNumber_BDSKY_Serial.7, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.8, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                  Delta_Re_sky$reproductiveNumber_BDSKY_Serial.10))
DeltaRe_10$ID <- factor(c('Delta'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/Delta/15_dimensions")
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

### dataframe with Re values and hpds 

Delta_times <- format(date_decimal(Delta_times), "%Y-%m-%d")

Delta_Re <- as.data.frame(t(Delta_Re_gridded_hpd))
Delta_Re <- cbind(Delta_Re, Delta_times)
names(Delta_Re) <-c("Lower Bound","mean","Upper Bound","date")

Delta_Re <- as.data.frame(t(Delta_Re_hpd))
Delta_Re$ID <- factor(c('Delta'))

DeltaRe_15 <- as.data.frame(cbind(Delta_Re_sky$reproductiveNumber_BDSKY_Serial.1, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.2, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                  Delta_Re_sky$reproductiveNumber_BDSKY_Serial.4, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.5, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                  Delta_Re_sky$reproductiveNumber_BDSKY_Serial.7, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.8, Delta_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                  Delta_Re_sky$reproductiveNumber_BDSKY_Serial.10))
DeltaRe_15$ID <- factor(c('Delta'))

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

### dataframe with Re values and hpds 

P1_times <- format(date_decimal(P1_times), "%Y-%m-%d")

P1_Re <- as.data.frame(t(P1_Re_gridded_hpd))
P1_Re <- cbind(P1_Re, P1_times)
names(P1_Re) <-c("Lower Bound","mean","Upper Bound","date")

P1_Re <- as.data.frame(t(P1_Re_hpd))
P1_Re$ID <- factor(c('Gamma'))

P1_Re_10 <- as.data.frame(cbind(P1_Re_sky$reproductiveNumber_BDSKY_Serial.1, P1_Re_sky$reproductiveNumber_BDSKY_Serial.2, P1_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
               P1_Re_sky$reproductiveNumber_BDSKY_Serial.4, P1_Re_sky$reproductiveNumber_BDSKY_Serial.5, P1_Re_sky$reproductiveNumber_BDSKY_Serial.6,
               P1_Re_sky$reproductiveNumber_BDSKY_Serial.7, P1_Re_sky$reproductiveNumber_BDSKY_Serial.8, P1_Re_sky$reproductiveNumber_BDSKY_Serial.9,
               P1_Re_sky$reproductiveNumber_BDSKY_Serial.10))
P1_Re_10$ID <- factor(c('Gamma'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/P_1/15_dimensions")
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

### dataframe with Re values and hpds 

P1_times <- format(date_decimal(P1_times), "%Y-%m-%d")

P1_Re <- as.data.frame(t(P1_Re_gridded_hpd))
P1_Re <- cbind(P1_Re, P1_times)
names(P1_Re) <-c("Lower Bound","mean","Upper Bound","date")

P1_Re <- as.data.frame(t(P1_Re_hpd))
P1_Re$ID <- factor(c('Gamma'))

P1_Re_15 <- as.data.frame(cbind(P1_Re_sky$reproductiveNumber_BDSKY_Serial.1, P1_Re_sky$reproductiveNumber_BDSKY_Serial.2, P1_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                P1_Re_sky$reproductiveNumber_BDSKY_Serial.4, P1_Re_sky$reproductiveNumber_BDSKY_Serial.5, P1_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                P1_Re_sky$reproductiveNumber_BDSKY_Serial.7, P1_Re_sky$reproductiveNumber_BDSKY_Serial.8, P1_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                P1_Re_sky$reproductiveNumber_BDSKY_Serial.10))
P1_Re_15$ID <- factor(c('Gamma'))

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

### dataframe with Re values and hpds 

B_1_1_7_times <- format(date_decimal(B_1_1_7_times), "%Y-%m-%d")

B_1_1_7_Re <- as.data.frame(t(B_1_1_7_Re_gridded_hpd))
B_1_1_7_Re <- cbind(B_1_1_7_Re, B_1_1_7_times)
names(B_1_1_7_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_1_7_Re <- as.data.frame(t(B_1_1_7_Re_hpd))
B_1_1_7_Re$ID <- factor(c('Alpha'))

B117_Re_10 <- as.data.frame(cbind(B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                 B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                 B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                 B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B117_Re_10$ID <- factor(c('Alpha'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_1_7/15_dimensions")
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

### dataframe with Re values and hpds 

B_1_1_7_times <- format(date_decimal(B_1_1_7_times), "%Y-%m-%d")

B_1_1_7_Re <- as.data.frame(t(B_1_1_7_Re_gridded_hpd))
B_1_1_7_Re <- cbind(B_1_1_7_Re, B_1_1_7_times)
names(B_1_1_7_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_1_7_Re <- as.data.frame(t(B_1_1_7_Re_hpd))
B_1_1_7_Re$ID <- factor(c('Alpha'))

B117_Re_15 <- as.data.frame(cbind(B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                  B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                  B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                  B_1_1_7_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B117_Re_15$ID <- factor(c('Alpha'))

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

### dataframe with Re values and hpds 

C37_times <- format(date_decimal(C37_times), "%Y-%m-%d")

C37_Re <- as.data.frame(t(C37_Re_gridded_hpd))
C37_Re <- cbind(C37_Re, C37_times)
names(C37_Re) <-c("Lower Bound","mean","Upper Bound","date")

C37_Re <- as.data.frame(t(C37_Re_hpd))
C37_Re$ID <- factor(c('Lambda'))

C37Re_10 <- as.data.frame(cbind(C37_Re_sky$reproductiveNumber_BDSKY_Serial.1, C37_Re_sky$reproductiveNumber_BDSKY_Serial.2, C37_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
               C37_Re_sky$reproductiveNumber_BDSKY_Serial.4, C37_Re_sky$reproductiveNumber_BDSKY_Serial.5, C37_Re_sky$reproductiveNumber_BDSKY_Serial.6,
               C37_Re_sky$reproductiveNumber_BDSKY_Serial.7, C37_Re_sky$reproductiveNumber_BDSKY_Serial.8, C37_Re_sky$reproductiveNumber_BDSKY_Serial.9,
               C37_Re_sky$reproductiveNumber_BDSKY_Serial.10))
C37Re_10$ID <- factor(c('Lambda'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/C_37/15_dimensions")
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

### dataframe with Re values and hpds 

C37_times <- format(date_decimal(C37_times), "%Y-%m-%d")

C37_Re <- as.data.frame(t(C37_Re_gridded_hpd))
C37_Re <- cbind(C37_Re, C37_times)
names(C37_Re) <-c("Lower Bound","mean","Upper Bound","date")

C37_Re <- as.data.frame(t(C37_Re_hpd))
C37_Re$ID <- factor(c('Lambda'))

C37Re_15 <- as.data.frame(cbind(C37_Re_sky$reproductiveNumber_BDSKY_Serial.1, C37_Re_sky$reproductiveNumber_BDSKY_Serial.2, C37_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                C37_Re_sky$reproductiveNumber_BDSKY_Serial.4, C37_Re_sky$reproductiveNumber_BDSKY_Serial.5, C37_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                C37_Re_sky$reproductiveNumber_BDSKY_Serial.7, C37_Re_sky$reproductiveNumber_BDSKY_Serial.8, C37_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                C37_Re_sky$reproductiveNumber_BDSKY_Serial.10))
C37Re_15$ID <- factor(c('Lambda'))

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

### dataframe with Re values and hpds 

Omi_times <- format(date_decimal(Omi_times), "%Y-%m-%d")

Omi_Re <- as.data.frame(t(Omi_Re_gridded_hpd))
Omi_Re <- cbind(Omi_Re, Omi_times)
names(Omi_Re) <-c("Lower Bound","mean","Upper Bound","date")

Omi_Re <- as.data.frame(t(Omi_Re_hpd))
Omi_Re$ID <- factor(c('Omicron'))


OmiRe_10 <- as.data.frame(cbind(Omi_Re_sky$reproductiveNumber_BDSKY_Serial.1, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.2, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
               Omi_Re_sky$reproductiveNumber_BDSKY_Serial.4, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.5, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.6,
               Omi_Re_sky$reproductiveNumber_BDSKY_Serial.7, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.8, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.9,
               Omi_Re_sky$reproductiveNumber_BDSKY_Serial.10))
OmiRe_10$ID <- factor(c('Omicron'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/Omicron/15_dimensions")
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

### dataframe with Re values and hpds 

Omi_times <- format(date_decimal(Omi_times), "%Y-%m-%d")

Omi_Re <- as.data.frame(t(Omi_Re_gridded_hpd))
Omi_Re <- cbind(Omi_Re, Omi_times)
names(Omi_Re) <-c("Lower Bound","mean","Upper Bound","date")

Omi_Re <- as.data.frame(t(Omi_Re_hpd))
Omi_Re$ID <- factor(c('Omicron'))


OmiRe_15 <- as.data.frame(cbind(Omi_Re_sky$reproductiveNumber_BDSKY_Serial.1, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.2, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                Omi_Re_sky$reproductiveNumber_BDSKY_Serial.4, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.5, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                Omi_Re_sky$reproductiveNumber_BDSKY_Serial.7, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.8, Omi_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                Omi_Re_sky$reproductiveNumber_BDSKY_Serial.10))
OmiRe_15$ID <- factor(c('Omicron'))

#### B_1_621

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_621/10_dimensions")
fname <- "BDSKL2_COMB.log"    
B_1_621_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_621_Re_sky    <- getSkylineSubset(B_1_621_lf, "reproductiveNumber")
B_1_621_Re_hpd    <- getMatrixHPD(B_1_621_Re_sky)
B_1_621_delta_hpd <- getHPD(B_1_621_lf$becomeUninfectiousRate)

plotSkyline(1:10, B_1_621_Re_hpd, type='step', ylab="R")


# plot the smooth skyline 2021-01-11 to 2021-12-06 
period = 0.95
timegrid       <- seq(0,period,length.out=47)
B_1_621_Re_gridded     <- gridSkyline(B_1_621_Re_sky, B_1_621_lf$origin, timegrid)
B_1_621_Re_gridded_hpd <- getMatrixHPD(B_1_621_Re_gridded)

# last sample date 2021-12-06 (329)/365
date = 329/365 + 2021
B_1_621_times <- date - timegrid


### dataframe with Re values and hpds 

B_1_621_times <- format(date_decimal(B_1_621_times), "%Y-%m-%d")

B_1_621_Re <- as.data.frame(t(B_1_621_Re_gridded_hpd))
B_1_621_Re <- cbind(B_1_621_Re, B_1_621_times)
names(B_1_621_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_621_Re <- as.data.frame(t(B_1_621_Re_hpd))
B_1_621_Re$ID <- factor(c('Mu'))

B_1_621_Re_10 <- as.data.frame(cbind(B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                             B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                             B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                             B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B_1_621_Re_10$ID <- factor(c('Mu'))

setwd("/Volumes/GoogleDrive-101426876281184082428/My Drive/Phd/SARS-CoV-2_Colombia_Phylodinamics-/Colombia/BDSK/B_1_621/15_dimensions")
fname <- "BDSKL2_COMB.log"    
B_1_621_lf    <- readLogfile(fname, burnin=0.1)

# extract the HPDs of Re and the becoming uninfectious rate:
B_1_621_Re_sky    <- getSkylineSubset(B_1_621_lf, "reproductiveNumber")
B_1_621_Re_hpd    <- getMatrixHPD(B_1_621_Re_sky)
B_1_621_delta_hpd <- getHPD(B_1_621_lf$becomeUninfectiousRate)

plotSkyline(1:10, B_1_621_Re_hpd, type='step', ylab="R")


# plot the smooth skyline 2021-01-11 to 2021-12-06 
period = 0.95
timegrid       <- seq(0,period,length.out=47)
B_1_621_Re_gridded     <- gridSkyline(B_1_621_Re_sky, B_1_621_lf$origin, timegrid)
B_1_621_Re_gridded_hpd <- getMatrixHPD(B_1_621_Re_gridded)

# last sample date 2021-12-06 (329)/365
date = 329/365 + 2021
B_1_621_times <- date - timegrid


### dataframe with Re values and hpds 

B_1_621_times <- format(date_decimal(B_1_621_times), "%Y-%m-%d")

B_1_621_Re <- as.data.frame(t(B_1_621_Re_gridded_hpd))
B_1_621_Re <- cbind(B_1_621_Re, B_1_621_times)
names(B_1_621_Re) <-c("Lower Bound","mean","Upper Bound","date")

B_1_621_Re <- as.data.frame(t(B_1_621_Re_hpd))
B_1_621_Re$ID <- factor(c('Mu'))

B_1_621_Re_15 <- as.data.frame(cbind(B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.1, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.2, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.3, 
                                     B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.4, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.5, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.6,
                                     B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.7, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.8, B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.9,
                                     B_1_621_Re_sky$reproductiveNumber_BDSKY_Serial.10))
B_1_621_Re_15$ID <- factor(c('Mu'))

#plotSkyline(B_1_621_times, B_1_621_Re_gridded_hpd, type='smooth', xlab="Time", ylab="R")


###############################

### Boxplot of all Re 


Re_10 <- rbind(B1Re_10, B111Re_10, B11348Re_10, B_1_420_Re, C37Re_10, B117_Re_10, B_1_621_Re_10, P1_Re_10, DeltaRe_10, OmiRe_10)

Re_15 <- rbind(B1Re_15, B111Re_15, B11348Re_15, B_1_420_Re, C37Re_15, B117_Re_15, B_1_621_Re_15, P1_Re_15, DeltaRe_15, OmiRe_10) 



dimensions_10 <- ggplot(Re_10, aes(x=ID, y=V2, color=ID)) + 
  geom_boxplot()+
  xlab("Variants") +
  ylab ("Re")+
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
  geom_hline(yintercept =1, colour = "red", lty=2, size=0.5, alpha=0.5)+
  theme_set(theme_classic(base_size=30))+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 25, angle=-90), axis.text.y = element_text(size = 25),
        axis.title.y=element_text(face="italic"),
        axis.title = element_text(size = 35), plot.title = element_text(face="bold", size=35))+
  ggtitle( "A.  10 dimensions")


dimensions_15 <- ggplot(Re_15, aes(x=ID, y=V2, color=ID)) + 
  geom_boxplot()+
  xlab("Variants") +
  ylab ("Re")+
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
  geom_hline(yintercept =1, colour = "red", lty=2, size=0.5, alpha=0.5)+
  theme_set(theme_classic(base_size=30))+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 25, angle=-90), axis.text.y = element_text(size = 25),
        axis.title.y=element_text(face="italic"),
        axis.title = element_text(size = 35), plot.title = element_text(face="bold", size=35))+
  ggtitle( "B.  14 dimensions")


  
  
pdf ('S5Fig.pdf',width=20,height=10)
  #par(mfrow=c(2,4))
  
  ggarrange(dimensions_10,dimensions_15, ncol = 2, nrow = 1)
  
dev.off()
  









