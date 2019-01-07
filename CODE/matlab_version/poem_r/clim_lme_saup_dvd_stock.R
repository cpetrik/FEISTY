################################################################################

# Comparisons of catch and frac pelagic 

################################################################################

rm(list=ls())

library(betareg)
library(visreg)
# library(ggplot2)
# library(gridExtra)
# library(arm)
# library(effects)
# library(plyr)
library(lattice) #multipanel plots
library(mgcv) #gam
# library(MuMIn)
# library(mnormt)
# library(modEvA)


#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
cfile = "Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100"
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/")

# load data
#all LMEs
AL <- read.csv(paste0("LME_DvD_SAU_vals_",cfile,".csv"),sep=",",header = T,stringsAsFactors = T)
#not LELC
NL <- read.csv(paste0("LME_DvD_SAU_notLELC_vals_",cfile,".csv"),sep=",",header = T,stringsAsFactors = T)

# NaNs
sid <- 64
did <- c(62,64,65,66)
did2 <- 45

### LINEAR MODEL ------------------------------------------------------------------
## FRACTION P VS. D
#All LMEs
Sall <- AL[-sid,10:12]
Dall <- AL[-did,10:12]

cor(Sall$PPD,Sall$SPD) #0.1298395
mPDs <- betareg(PPD ~ SPD -1, data=Sall)
summary(mPDs)
visreg(mPDs)

cor(Dall$PPD,Dall$DPD) #0.3086556
mPDd <- betareg(PPD ~ DPD -1, data=Dall)
summary(mPDd)
visreg(mPDd)

#notLELC
Dall2 <- NL[-did2,10:12]

cor(NL$PPD,NL$SPD) #0.2237171
lPDs <- betareg(PPD ~ SPD -1, data=NL)
summary(lPDs)
visreg(lPDs)

cor(Dall2$PPD,Dall2$DPD) #0.2713948
lPDd <- betareg(PPD ~ DPD -1, data=Dall2)
summary(lPDd)
visreg(lPDd)


