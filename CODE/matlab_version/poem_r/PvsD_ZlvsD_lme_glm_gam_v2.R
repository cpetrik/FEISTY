################################################################################

# GLMs of frac pelagic as a function of 
# Zl:D ratio, temperature, NPP, coastal area
# New parameters kt=0.0805, BE=0.075

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
cfile = "Dc_enc70-b200_m4-b175-k08_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100"
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/")

# load data
ZB <- read.csv("LME_ZBratios_clim_fished_All_fish03.csv",sep=",",header = T,stringsAsFactors = T)
ZB$logRatZDet <- log10(ZB$RatZDet)
ZB$logRatZB <- log10(ZB$RatZB)
ZB$logRatZlDet <- log10(ZB$RatZlDet)
ZB$logRatZlB <- log10(ZB$RatZlB)
ZB$logNPP <- log10(ZB$NPP/365)

### DATA EXPLORATION ------------------------------------------------------------------
# Outliers
dat <- c("LME_ptemp","FracPD","FracPF","FracLM","LME_depth","LME_Frac200",
         "logRatZDet","logRatZB","logRatZlDet","logRatZlB","logNPP")
Mydotplot(as.matrix(ZB[,dat]))
#npp has a few low values

# Collinearity
pairs(ZB[,dat],lower.panel = panel.cor)
#lme depth >=0.7 corr to ratios
#npp 0.5 with ZlDet and 0.7 with ZlB

# Linear relationships?
vars <- c("LME_ptemp","LME_Frac200","logRatZDet","logRatZB","logRatZlDet",
          "logRatZlB","logNPP")
Myxyplot(ZB,vars,"FracPD") #temp, ZDet curves
Myxyplot(ZB,vars,"FracPF") #temp, ZDet curves
Myxyplot(ZB,vars,"FracLM") #all or none?

### LINEAR MODEL ------------------------------------------------------------------
## FRACTION P VS. D
mPDf <- betareg(FracPD ~ LME_ptemp + logRatZlDet + LME_Frac200, data=ZB)
summary(mPDf)
mPDd <- betareg(FracPD ~ LME_ptemp + logRatZlDet + LME_depth, data=ZB)
summary(mPDd)
mPDt <- betareg(FracPD ~ LME_ptemp, data=ZB)
summary(mPDt) 
mPDr <- betareg(FracPD ~ logRatZlDet, data=ZB)
summary(mPDr)
mPDf2 <- betareg(FracPD ~ LME_Frac200, data=ZB)
summary(mPDf2)
mPDd2 <- betareg(FracPD ~ LME_depth, data=ZB) #depth better than frac200
summary(mPDd2)
mPDn <- betareg(FracPD ~ logNPP, data=ZB)
summary(mPDn)

par(mfrow=c(2,2))
visreg(mPDf)
par(mfrow=c(2,2))
visreg(mPDd)
par(mfrow=c(1,1))
visreg(mPDr)


## FRACTION P VS. F
mPFf <- betareg(FracPF ~ LME_ptemp + logRatZlDet + LME_Frac200, data=ZB)
summary(mPFf)
mPFd <- betareg(FracPF ~ LME_ptemp + logRatZlDet + LME_depth, data=ZB)
summary(mPFd)
mPFt <- betareg(FracPF ~ LME_ptemp, data=ZB)
summary(mPFt) 
mPFr <- betareg(FracPF ~ logRatZlDet, data=ZB)
summary(mPFr)
mPFf2 <- betareg(FracPF ~ LME_Frac200, data=ZB) #frac better than depth
summary(mPFf2)
mPFd2 <- betareg(FracPF ~ LME_depth, data=ZB)
summary(mPFd2)
mPFn <- betareg(FracPF ~ logNPP, data=ZB)
summary(mPFn)

par(mfrow=c(2,2))
visreg(mPFf)
par(mfrow=c(2,2))
visreg(mPFd)
par(mfrow=c(1,1))
visreg(mPFr)


## FRACTION L VS. M
mLMf <- betareg(FracLM ~ LME_ptemp + logRatZlDet + LME_Frac200, data=ZB)
summary(mLMf)
mLMd <- betareg(FracLM ~ LME_ptemp + logRatZlDet + LME_depth, data=ZB)
summary(mLMd)
mLMt <- betareg(FracLM ~ LME_ptemp, data=ZB)
summary(mLMt) 
mLMr <- betareg(FracLM ~ logRatZlDet, data=ZB)
summary(mLMr)
mLMf2 <- betareg(FracLM ~ LME_Frac200, data=ZB) 
summary(mLMf2)
mLMd2 <- betareg(FracLM ~ LME_depth, data=ZB) #depth better than frac
summary(mLMd2)
mLMn <- betareg(FracLM ~ logNPP, data=ZB)
summary(mLMn)

par(mfrow=c(2,2))
visreg(mLMf)
par(mfrow=c(2,2))
visreg(mLMd)
par(mfrow=c(1,1))
visreg(mLMr)


## Model validation using ZlDet
# PD
EPD <- resid(mPDf)
FPD <- fitted(mPDf)
pdf(paste0(fpath,"LME_PDfrac_ZlDet_ptemp_frac200_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#looks like some nonlinearity in the residuals and definitely outliers

EPD <- resid(mPDr)
FPD <- fitted(mPDr)
pdf(paste0(fpath,"LME_PDfrac_ZlDet_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#no nonlinearity, just outliers

# PF
EPF <- resid(mPFf)
FPF <- fitted(mPFf)
pdf(paste0(fpath,"LME_PFfrac_ZlDet_ptemp_frac200_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
# Looks like PD

EPF <- resid(mPFr)
FPF <- fitted(mPFr)
pdf(paste0(fpath,"LME_PFfrac_ZlDet_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#okay

# LM
ELM <- resid(mLMf)
FLM <- fitted(mLMf)
pdf(paste0(fpath,"LME_LMfrac_ZlDet_ptemp_frac200_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#def nonlinear

ELM <- resid(mLMr)
FLM <- fitted(mLMr)
pdf(paste0(fpath,"LME_LMfrac_ZlDet_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have neg relat with temp



### ADDITIVE MODEL ------------------------------------------------------------------
## NEED TO GIVE IT BETA DISTRIBUTION (family=beta)?

## Model selection using ZlDet
gPD <- gam(FracPD ~ s(LME_ptemp,k=4) + s(logRatZlDet,k=4) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPD) #86.8% deviance
gPDl <- gam(FracPD ~ s(LME_ptemp,k=4) + (logRatZlDet) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPDl) #84.9% dev
gPDr <- gam(FracPD ~ s(logRatZlDet,k=4), data=ZB, family = betar)
summary(gPDr) #69.2%
gPDt <- gam(FracPD ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gPDt) #66.2%
gPDf <- gam(FracPD ~ s(LME_Frac200,k=4), data=ZB, family = betar) #edf=1 linear
summary(gPDf) #26.6%
gPDn <- gam(FracPD ~ s(logNPP,k=4), data=ZB, family = betar) 
summary(gPDn) #71%
AIC(gPD,gPDr,gPDl,gPDt,gPDf,gPDn) #gPD > gPDl > gPDn > gPDr
AIC(gPD,mPDf) #gPD

gPF <- gam(FracPF ~ s(LME_ptemp,k=4) + s(logRatZlDet,k=4) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPF) #71.1%
gPFl <- gam(FracPF ~ s(LME_ptemp,k=4) + (logRatZlDet) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPFl) #68%
gPFr <- gam(FracPF ~ s(logRatZlDet,k=4), data=ZB, family = betar) #edf=1 linear
summary(gPFr) #46.9% dev
gPFt <- gam(FracPF ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gPFt) #43.7%
gPFf <- gam(FracPF ~ s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPFf) #10.4%
gPFn <- gam(FracPF ~ s(logNPP,k=4), data=ZB, family = betar) 
summary(gPFn) #58.6%
AIC(gPF,gPFr,gPFl,gPFt,gPFf,gPFn) #gPF > gPFl > gPFn > gPFr
AIC(gPF,mPFf) #gPF

gLM <- gam(FracLM ~ s(LME_ptemp,k=4) + s(logRatZlDet,k=4) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gLM) #70.8%
gLMl <- gam(FracLM ~ s(LME_ptemp,k=4) + (logRatZlDet) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gLMl) #70.8%
gLMr <- gam(FracLM ~ s(logRatZlDet,k=4), data=ZB, family = betar)
summary(gLMr) #23.8% dev
gLMt <- gam(FracLM ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gLMt) #57.3%
gLMf <- gam(FracLM ~ s(LME_Frac200,k=4), data=ZB, family = betar) #edf=1 linear
summary(gLMf) #12.7%
gLMn <- gam(FracLM ~ s(logNPP,k=4), data=ZB, family = betar) 
summary(gLMn) #22.8%
AIC(gLM,gLMr,gLMl,gLMt,gLMf,gLMn) #gLM=gLMl > gLMt > gLMr 
AIC(gLM,mLMf) #gLM

par(mfrow=c(2,2))
visreg(gPD)
par(mfrow=c(2,2))
visreg(gPF)
par(mfrow=c(2,2))
visreg(gLM)
#All have a kind of quadratic, concave-down relat with temp

par(mfrow=c(2,2))
visreg(gPDr)
visreg(gPFr)
visreg(gLMr)

#Zl:D rat linear
par(mfrow=c(1,3))
visreg(gPDl)
par(mfrow=c(1,3))
visreg(gPFl)
par(mfrow=c(1,3))
visreg(gLMl)
#All have a kind of quadratic, concave-down relat with temp


## Model validation using ZlDet
# PD
EPD <- resid(gPD)
FPD <- fitted(gPD)
pdf(paste0(fpath,"LME_PDfrac_ZlDet_ptemp_frac200_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#no nonlinearity, just outliers

EPD <- resid(gPDr)
FPD <- fitted(gPDr)
pdf(paste0(fpath,"LME_PDfrac_ZlDet_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#no nonlinearity, just outliers

# PF
EPF <- resid(gPF)
FPF <- fitted(gPF)
pdf(paste0(fpath,"LME_PFfrac_ZlDet_ptemp_frac200_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#OK

EPF <- resid(gPFr)
FPF <- fitted(gPFr)
pdf(paste0(fpath,"LME_PFfrac_ZlDet_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have neg relat with temp

# LM
ELM <- resid(gLM)
FLM <- fitted(gLM)
pdf(paste0(fpath,"LME_LMfrac_ZlDet_ptemp_frac200_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#OK

ELM <- resid(gLMr)
FLM <- fitted(gLMr)
pdf(paste0(fpath,"LME_LMfrac_ZlDet_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have strong relat with temp

##Test knots
gPDr3 <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar)
gPDr5 <- gam(FracPD ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gPDr6 <- gam(FracPD ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gPDr3,gPDr,gPDr5,gPDr6) #3 by 0.02

gPFr3 <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar)
gPFr5 <- gam(FracPF ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gPFr6 <- gam(FracPF ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gPFr3,gPFr,gPFr5,gPFr6) #3 by 0.8

gLMr3 <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar)
gLMr5 <- gam(FracLM ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gLMr6 <- gam(FracLM ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gLMr3,gLMr,gLMr5,gLMr6) #practically identical

##Test link fun
gPDrP <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gPDrL <- gam(FracPD ~ s(logRatZlDet,k=5), data=ZB, family = betar(link="cloglog"))
gPDrC <- gam(FracPD ~ s(logRatZlDet,k=6), data=ZB, family = betar(link="cauchit"))
AIC(gPDrP,gPDr,gPDrL,gPDrC) #probit>logit by 0.01

gPFrP <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gPFrL <- gam(FracPF ~ s(logRatZlDet,k=5), data=ZB, family = betar(link="cloglog"))
gPFrC <- gam(FracPF ~ s(logRatZlDet,k=6), data=ZB, family = betar(link="cauchit"))
AIC(gPFrP,gPFr,gPFrL,gPFrC) #probit>logit by 0.4

gLMrP <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gLMrL <- gam(FracLM ~ s(logRatZlDet,k=5), data=ZB, family = betar(link="cloglog"))
gLMrC <- gam(FracLM ~ s(logRatZlDet,k=6), data=ZB, family = betar(link="cauchit"))
AIC(gLMrP,gLMr,gLMrL,gLMrC) #cloglog by 0.006


### UPDATE TO 3 KNOTS --------------------------------------------------
## Model selection using ZlDet
g3PD <- gam(FracPD ~ s(LME_ptemp,k=3) + s(logRatZlDet,k=3) + s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3PD) #83.7% deviance
g3PDl <- gam(FracPD ~ s(LME_ptemp,k=3) + (logRatZlDet) + s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3PDl) #83.7% dev
g3PDr <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar)
summary(g3PDr) #68.6%
g3PDt <- gam(FracPD ~ s(LME_ptemp,k=3), data=ZB, family = betar)
summary(g3PDt) #65.8%
g3PDf <- gam(FracPD ~ s(LME_Frac200,k=3), data=ZB, family = betar) #edf=1 linear
summary(g3PDf) #26.6%
g3PDn <- gam(FracPD ~ s(logNPP,k=3), data=ZB, family = betar) 
summary(g3PDn) #69.9%
AIC(g3PD,g3PDr,g3PDl,g3PDt,g3PDf,g3PDn) #g3PD=g3PDl > g3PDn > g3PDr
AIC(g3PD,mPDf) #g3PD

g3PF <- gam(FracPF ~ s(LME_ptemp,k=3) + s(logRatZlDet,k=3) + s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3PF) #64.6%
g3PFl <- gam(FracPF ~ s(LME_ptemp,k=3) + (logRatZlDet) + s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3PFl) #64.6%
g3PFr <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar) #edf=1 linear
summary(g3PFr) #46.8% dev
g3PFt <- gam(FracPF ~ s(LME_ptemp,k=3), data=ZB, family = betar)
summary(g3PFt) #42.1%
g3PFf <- gam(FracPF ~ s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3PFf) #10.4%
g3PFn <- gam(FracPF ~ s(logNPP,k=3), data=ZB, family = betar) 
summary(g3PFn) #58.5%
AIC(g3PF,g3PFr,g3PFl,g3PFt,g3PFf,g3PFn) #g3PF > g3PFl > g3PFn > g3PFr
AIC(g3PF,mPFf) #g3PF

g3LM <- gam(FracLM ~ s(LME_ptemp,k=3) + s(logRatZlDet,k=3) + s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3LM) #68.7%
g3LMl <- gam(FracLM ~ s(LME_ptemp,k=3) + (logRatZlDet) + s(LME_Frac200,k=3), data=ZB, family = betar)
summary(g3LMl) #68.7%
g3LMr <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar)
summary(g3LMr) #23.8% dev
g3LMt <- gam(FracLM ~ s(LME_ptemp,k=3), data=ZB, family = betar)
summary(g3LMt) #55.6%
g3LMf <- gam(FracLM ~ s(LME_Frac200,k=3), data=ZB, family = betar) #edf=1 linear
summary(g3LMf) #12.7%
g3LMn <- gam(FracLM ~ s(logNPP,k=3), data=ZB, family = betar) 
summary(g3LMn) #4.49%
AIC(g3LM,g3LMr,g3LMl,g3LMt,g3LMf,g3LMn) #g3LM=g3LMl > g3LMt > g3LMr 
AIC(g3LM,mLMf) #g3LM

par(mfrow=c(2,2))
visreg(g3PD)
par(mfrow=c(2,2))
visreg(g3PF)
par(mfrow=c(2,2))
visreg(g3LM)
#All have a kind of quadratic, concave-down relat with temp

par(mfrow=c(2,2))
visreg(g3PDr)
visreg(g3PFr)
visreg(g3LMr)

#Zl:D rat linear
par(mfrow=c(1,3))
visreg(g3PDl)
par(mfrow=c(1,3))
visreg(g3PFl)
par(mfrow=c(1,3))
visreg(g3LMl)
#All have a kind of quadratic, concave-down relat with temp




### VISREG FIGURES --------------------------------------------------
#PD
pdf(paste0(fpath,"LME_PDfrac_ZlDet_glm_gam.pdf"))
par(mfrow=c(2,3))
visreg(mPDf,ylab="frac P vs. D")
visreg(gPD,ylab="frac P vs. D",scale="response")
dev.off()

#PF
pdf(paste0(fpath,"LME_PFfrac_ZlDet_glm_gam.pdf"))
par(mfrow=c(2,3))
visreg(mPFf,ylab="frac P vs. F")
visreg(gPF,ylab="frac P vs. F",scale="response")
dev.off()

#LM
pdf(paste0(fpath,"LME_LMfrac_ZlDet_glm_gam.pdf"))
par(mfrow=c(2,3))
visreg(mLMf,ylab="frac L vs. M")
visreg(gLM,ylab="frac L vs. M",scale="response")
dev.off()

#Just Zl:Det GLM
pdf(paste0(fpath,"LME_All_ZlDet_glm.pdf"))
par(mfrow=c(2,2))
visreg(mPDr,ylab="frac P vs. D")
visreg(mPFr,ylab="frac P vs. F")
visreg(mLMr,ylab="frac L vs. M")
dev.off()

#Just Zl:Det GAM
pdf(paste0(fpath,"LME_All_ZlDet_gam.pdf"))
par(mfrow=c(2,2))
visreg(gPDr,ylab="frac P vs. D",scale="response")
visreg(gPFr,ylab="frac P vs. F",scale="response")
visreg(gLMr,ylab="frac L vs. M",scale="response")
dev.off()

#Just pTemp GLM
pdf(paste0(fpath,"LME_All_ptemp_glm.pdf"))
par(mfrow=c(2,2))
visreg(mPDt,ylab="frac P vs. D")
visreg(mPFt,ylab="frac P vs. F")
visreg(mLMt,ylab="frac L vs. M")
dev.off()

#Just pTemp GAM
pdf(paste0(fpath,"LME_All_ptemp_gam.pdf"))
par(mfrow=c(2,2))
visreg(gPDt,ylab="frac P vs. D",scale="response")
visreg(gPFt,ylab="frac P vs. F",scale="response")
visreg(gLMt,ylab="frac L vs. M",scale="response")
dev.off()

#Just frac200 GLM
pdf(paste0(fpath,"LME_All_frac200_glm.pdf"))
par(mfrow=c(2,2))
visreg(mPDf2,ylab="frac P vs. D")
visreg(mPFf2,ylab="frac P vs. F")
visreg(mLMf2,ylab="frac L vs. M")
dev.off()

#Just frac200 GAM
pdf(paste0(fpath,"LME_All_frac200_gam.pdf"))
par(mfrow=c(2,2))
visreg(gPDf,ylab="frac P vs. D",scale="response")
visreg(gPFf,ylab="frac P vs. F",scale="response")
visreg(gLMf,ylab="frac L vs. M",scale="response")
dev.off()

#Just NPP GLM
pdf(paste0(fpath,"LME_All_npp_glm.pdf"))
par(mfrow=c(2,2))
visreg(mPDn,ylab="frac P vs. D")
visreg(mPFn,ylab="frac P vs. F")
visreg(mLMn,ylab="frac L vs. M")
dev.off()

#Just NPP GAM
pdf(paste0(fpath,"LME_All_npp_gam.pdf"))
par(mfrow=c(2,2))
visreg(gPDn,ylab="frac P vs. D",scale="response")
visreg(gPFn,ylab="frac P vs. F",scale="response")
visreg(gLMn,ylab="frac L vs. M",scale="response")
dev.off()



### PREDICT VALUES --------------------------------------------------
#Z:Det
logRatZlDet <- data.frame(logRatZlDet = seq(from=min(ZB$logRatZlDet),to=max(ZB$logRatZlDet),length=100))
PD <- predict(gPDr, newdata=logRatZlDet, type = "response", se=TRUE)
PF <- predict(gPFr, newdata=logRatZlDet, type = "response", se=TRUE)
LM <- predict(gLMr, newdata=logRatZlDet, type = "response", se=TRUE)

Dfit <- as.data.frame(logRatZlDet)
Dfit[,2] <- as.data.frame(PD$fit)
Dfit[,3] <- as.data.frame(PD$se.fit)
Dfit[,4] <- as.data.frame(PF$fit)
Dfit[,5] <- as.data.frame(PF$se.fit)
Dfit[,6] <- as.data.frame(LM$fit)
Dfit[,7] <- as.data.frame(LM$se.fit)
names(Dfit) <- c("logZlDet","PDfit","PDse","PFfit","PFse","LMfit","LMse")

write.table(Dfit,"ZlDet_gam_fit.csv",sep=",",row.names=F)


#Temp
LME_ptemp <- data.frame(LME_ptemp = seq(from=min(ZB$LME_ptemp),to=max(ZB$LME_ptemp),length=100))
tPD <- predict(gPDt, newdata=LME_ptemp, type = "response", se=TRUE)
tPF <- predict(gPFt, newdata=LME_ptemp, type = "response", se=TRUE)
tLM <- predict(gLMt, newdata=LME_ptemp, type = "response", se=TRUE)

Tfit <- as.data.frame(LME_ptemp)
Tfit[,2] <- as.data.frame(tPD$fit)
Tfit[,3] <- as.data.frame(tPD$se.fit)
Tfit[,4] <- as.data.frame(tPF$fit)
Tfit[,5] <- as.data.frame(tPF$se.fit)
Tfit[,6] <- as.data.frame(tLM$fit)
Tfit[,7] <- as.data.frame(tLM$se.fit)
names(Tfit) <- c("LME_ptemp","PDfit","PDse","PFfit","PFse","LMfit","LMse")

write.table(Tfit,"ptemp_gam_fit.csv",sep=",",row.names=F)


#Depth
LME_Frac200 <- data.frame(LME_Frac200 = seq(from=min(ZB$LME_Frac200),to=max(ZB$LME_Frac200),length=100))
fPD <- predict(gPDf, newdata=LME_Frac200, type = "response", se=TRUE)
fPF <- predict(gPFf, newdata=LME_Frac200, type = "response", se=TRUE)
fLM <- predict(gLMf, newdata=LME_Frac200, type = "response", se=TRUE)

Ffit <- as.data.frame(LME_Frac200)
Ffit[,2] <- as.data.frame(fPD$fit)
Ffit[,3] <- as.data.frame(fPD$se.fit)
Ffit[,4] <- as.data.frame(fPF$fit)
Ffit[,5] <- as.data.frame(fPF$se.fit)
Ffit[,6] <- as.data.frame(fLM$fit)
Ffit[,7] <- as.data.frame(fLM$se.fit)
names(Ffit) <- c("LME_Frac200","PDfit","PDse","PFfit","PFse","LMfit","LMse")

write.table(Ffit,"Frac200_gam_fit.csv",sep=",",row.names=F)


#NPP
logNPP <- data.frame(logNPP = seq(from=min(ZB$logNPP),to=max(ZB$logNPP),length=100))
nPD <- predict(gPDn, newdata=logNPP, type = "response", se=TRUE)
nPF <- predict(gPFn, newdata=logNPP, type = "response", se=TRUE)
nLM <- predict(gLMn, newdata=logNPP, type = "response", se=TRUE)

Nfit <- as.data.frame(logNPP)
Nfit[,2] <- as.data.frame(nPD$fit)
Nfit[,3] <- as.data.frame(nPD$se.fit)
Nfit[,4] <- as.data.frame(nPF$fit)
Nfit[,5] <- as.data.frame(nPF$se.fit)
Nfit[,6] <- as.data.frame(nLM$fit)
Nfit[,7] <- as.data.frame(nLM$se.fit)
names(Nfit) <- c("logNPP","PDfit","PDse","PFfit","PFse","LMfit","LMse")

write.table(Nfit,"npp_gam_fit.csv",sep=",",row.names=F)


visreg(gPDr)
visreg(gPDn)
visreg(gPDr,scale="response")
visreg(gPDn,scale="response")

### TABLE --------------------------------------------------
library("texreg")
screenreg(list(gPDr,gPDt,gPDf,gPDn))
screenreg(list(gPFr,gPFt,gPFf,gPFn))
screenreg(list(gLMr,gLMt,gLMf,gLMn))

htmlreg(list(gPDr,gPDt,gPDf,gPDn), file = "PDgam_table.doc", inline.css = FALSE, 
        doctype = TRUE, html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)
htmlreg(list(gPDr,gPDt,gPDf,gPDn), file = "PDgam_table.html")
htmlreg(list(gPFr,gPFt,gPFf,gPFn), file = "PFgam_table.html")
htmlreg(list(gLMr,gLMt,gLMf,gLMn), file = "LMgam_table.html")

