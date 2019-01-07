################################################################################

# GLMs of frac pelagic as a function of 
# Zl:D ratio, temperature, NPP, coastal area

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
Myxyplot(ZB,vars,"FracLM") #all ratios?

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
#resids have parabol relat with temp



### ADDITIVE MODEL ------------------------------------------------------------------
## NEED TO GIVE IT BETA DISTRIBUTION (family=beta)?

## Model selection using ZlDet
gPD <- gam(FracPD ~ s(LME_ptemp,k=4) + s(logRatZlDet,k=4) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPD) #81.4% deviance
gPDl <- gam(FracPD ~ s(LME_ptemp,k=4) + (logRatZlDet) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPDl)
gPDr <- gam(FracPD ~ s(logRatZlDet,k=4), data=ZB, family = betar)
summary(gPDr) #59.6% dev
gPDt <- gam(FracPD ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gPDt) #57.6%
gPDf <- gam(FracPD ~ s(LME_Frac200,k=4), data=ZB, family = betar) #edf=1 linear
summary(gPDf) #34.5%
gPDn <- gam(FracPD ~ s(logNPP,k=4), data=ZB, family = betar) 
summary(gPDn) #33.4%
AIC(gPD,gPDr,gPDl,gPDt,gPDf,gPDn) #gPD
AIC(gPD,mPDf) #gPD

gPF <- gam(FracPF ~ s(LME_ptemp,k=4) + s(logRatZlDet,k=4) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPF) #75%
gPFl <- gam(FracPF ~ s(LME_ptemp,k=4) + (logRatZlDet) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPFl)
gPFr <- gam(FracPF ~ s(logRatZlDet,k=4), data=ZB, family = betar) #edf=1 linear
summary(gPFr) #27.9% dev
gPFt <- gam(FracPF ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gPFt) #54.2%
gPFf <- gam(FracPF ~ s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gPFf) #23.4%
gPFn <- gam(FracPF ~ s(logNPP,k=4), data=ZB, family = betar) 
summary(gPFn) #14.5%
AIC(gPF,gPFr,gPFl,gPFt,gPFf,gPFn) #gPF
AIC(gPF,mPFf) #gPF

gLM <- gam(FracLM ~ s(LME_ptemp,k=4) + s(logRatZlDet,k=4) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gLM) #78.4%
gLMl <- gam(FracLM ~ s(LME_ptemp,k=4) + (logRatZlDet) + s(LME_Frac200,k=4), data=ZB, family = betar)
summary(gLMl)
gLMr <- gam(FracLM ~ s(logRatZlDet,k=4), data=ZB, family = betar)
summary(gLMr) #20.7% dev
gLMt <- gam(FracLM ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gLMt) #68.8%
gLMf <- gam(FracLM ~ s(LME_Frac200,k=4), data=ZB, family = betar) #edf=1 linear
summary(gLMf) #11.5%
gLMn <- gam(FracLM ~ s(logNPP,k=4), data=ZB, family = betar) 
summary(gLMn) #18.9%
AIC(gLM,gLMr,gLMl,gLMt,gLMf,gLMn) #gLM=gLMl
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
AIC(gPDr3,gPDr,gPDr5,gPDr6) #3 by 0.6

gPFr3 <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar)
gPFr5 <- gam(FracPF ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gPFr6 <- gam(FracPF ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gPFr3,gPFr,gPFr5,gPFr6) #3 by 0.2

gLMr3 <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar)
gLMr5 <- gam(FracLM ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gLMr6 <- gam(FracLM ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gLMr3,gLMr,gLMr5,gLMr6) #practically identical

##Test link fun
gPDrP <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gPDrL <- gam(FracPD ~ s(logRatZlDet,k=5), data=ZB, family = betar(link="cloglog"))
gPDrC <- gam(FracPD ~ s(logRatZlDet,k=6), data=ZB, family = betar(link="cauchit"))
AIC(gPDrP,gPDr,gPDrL,gPDrC) #logit by 0.9

gPFrP <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gPFrL <- gam(FracPF ~ s(logRatZlDet,k=5), data=ZB, family = betar(link="cloglog"))
gPFrC <- gam(FracPF ~ s(logRatZlDet,k=6), data=ZB, family = betar(link="cauchit"))
AIC(gPFrP,gPFr,gPFrL,gPFrC) #probit>logit by 0.2

gLMrP <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gLMrL <- gam(FracLM ~ s(logRatZlDet,k=5), data=ZB, family = betar(link="cloglog"))
gLMrC <- gam(FracLM ~ s(logRatZlDet,k=6), data=ZB, family = betar(link="cauchit"))
AIC(gLMrP,gLMr,gLMrL,gLMrC) #all equiv



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



### PREDICT VALUES --------------------------------------------------
#Z:Det
logRatZlDet <- data.frame(logRatZlDet = seq(from=min(ZB$logRatZlDet),to=max(ZB$logRatZlDet),length=100))
PD <- predict(gPDr, newdata=logRatZlDet, type = "response", se=TRUE)
PF <- predict(gPFr, newdata=logRatZlDet, type = "response", se=TRUE)
LM <- predict(gLMr, newdata=logRatZlDet, type = "response", se=TRUE)

PDfit <- as.data.frame(PD$fit)
PDfit[,2] <- as.data.frame(PD$se.fit)
PDfit[,3] <- as.data.frame(logRatZlDet)
names(PDfit) <- c("fit","se","logZlDet")

PFfit <- as.data.frame(PF$fit)
PFfit[,2] <- as.data.frame(PF$se.fit)
PFfit[,3] <- as.data.frame(logRatZlDet)
names(PFfit) <- c("fit","se","logZlDet")

LMfit <- as.data.frame(LM$fit)
LMfit[,2] <- as.data.frame(LM$se.fit)
LMfit[,3] <- as.data.frame(logRatZlDet)
names(LMfit) <- c("fit","se","logZlDet")

write.table(PDfit,"PD_ZlDet_gam_fit.csv",sep=",",row.names=F)
write.table(PFfit,"PF_ZlDet_gam_fit.csv",sep=",",row.names=F)
write.table(LMfit,"LM_ZlDet_gam_fit.csv",sep=",",row.names=F)


#Temp
LME_ptemp <- data.frame(LME_ptemp = seq(from=min(ZB$LME_ptemp),to=max(ZB$LME_ptemp),length=100))
tPD <- predict(gPDt, newdata=LME_ptemp, type = "response", se=TRUE)
tPF <- predict(gPFt, newdata=LME_ptemp, type = "response", se=TRUE)
tLM <- predict(gLMt, newdata=LME_ptemp, type = "response", se=TRUE)

tPDfit <- as.data.frame(tPD$fit)
tPDfit[,2] <- as.data.frame(tPD$se.fit)
tPDfit[,3] <- as.data.frame(LME_ptemp)
names(tPDfit) <- c("fit","se","LME_ptemp")

tPFfit <- as.data.frame(tPF$fit)
tPFfit[,2] <- as.data.frame(tPF$se.fit)
tPFfit[,3] <- as.data.frame(LME_ptemp)
names(tPFfit) <- c("fit","se","LME_ptemp")

tLMfit <- as.data.frame(tLM$fit)
tLMfit[,2] <- as.data.frame(tLM$se.fit)
tLMfit[,3] <- as.data.frame(LME_ptemp)
names(tLMfit) <- c("fit","se","LME_ptemp")

write.table(tPDfit,"PD_ptemp_gam_fit.csv",sep=",",row.names=F)
write.table(tPFfit,"PF_ptemp_gam_fit.csv",sep=",",row.names=F)
write.table(tLMfit,"LM_ptemp_gam_fit.csv",sep=",",row.names=F)


#Depth
LME_Frac200 <- data.frame(LME_Frac200 = seq(from=min(ZB$LME_Frac200),to=max(ZB$LME_Frac200),length=100))
fPD <- predict(gPDf, newdata=LME_Frac200, type = "response", se=TRUE)
fPF <- predict(gPFf, newdata=LME_Frac200, type = "response", se=TRUE)
fLM <- predict(gLMf, newdata=LME_Frac200, type = "response", se=TRUE)

fPDfit <- as.data.frame(fPD$fit)
fPDfit[,2] <- as.data.frame(fPD$se.fit)
fPDfit[,3] <- as.data.frame(LME_Frac200)
names(fPDfit) <- c("fit","se","LME_Frac200")

fPFfit <- as.data.frame(fPF$fit)
fPFfit[,2] <- as.data.frame(fPF$se.fit)
fPFfit[,3] <- as.data.frame(LME_Frac200)
names(fPFfit) <- c("fit","se","LME_Frac200")

fLMfit <- as.data.frame(fLM$fit)
fLMfit[,2] <- as.data.frame(fLM$se.fit)
fLMfit[,3] <- as.data.frame(LME_Frac200)
names(fLMfit) <- c("fit","se","LME_Frac200")

write.table(fPDfit,"PD_Frac200_gam_fit.csv",sep=",",row.names=F)
write.table(fPFfit,"PF_Frac200_gam_fit.csv",sep=",",row.names=F)
write.table(fLMfit,"LM_Frac200_gam_fit.csv",sep=",",row.names=F)


#NPP
logNPP <- data.frame(logNPP = seq(from=min(ZB$logNPP),to=max(ZB$logNPP),length=100))
nPD <- predict(gPDn, newdata=logNPP, type = "response", se=TRUE)
nPF <- predict(gPFn, newdata=logNPP, type = "response", se=TRUE)
nLM <- predict(gLMn, newdata=logNPP, type = "response", se=TRUE)

nPDfit <- as.data.frame(nPD$fit)
nPDfit[,2] <- as.data.frame(nPD$se.fit)
nPDfit[,3] <- as.data.frame(logNPP)
names(nPDfit) <- c("fit","se","logNPP")

nPFfit <- as.data.frame(nPF$fit)
nPFfit[,2] <- as.data.frame(nPF$se.fit)
nPFfit[,3] <- as.data.frame(logNPP)
names(nPFfit) <- c("fit","se","logNPP")

nLMfit <- as.data.frame(nLM$fit)
nLMfit[,2] <- as.data.frame(nLM$se.fit)
nLMfit[,3] <- as.data.frame(logNPP)
names(nLMfit) <- c("fit","se","logNPP")

write.table(nPDfit,"PD_npp_gam_fit.csv",sep=",",row.names=F)
write.table(nPFfit,"PF_npp_gam_fit.csv",sep=",",row.names=F)
write.table(nLMfit,"LM_npp_gam_fit.csv",sep=",",row.names=F)



### MANUAL FIGURES --------------------------------------------------
#PD
pdf(paste0(fpath,"LME_PDfrac_ZlDet_glm_gam_v2.pdf"))
par(mfrow=c(2,3))

dev.off()

#PF
pdf(paste0(fpath,"LME_PFfrac_ZlDet_glm_gam_v2.pdf"))
par(mfrow=c(2,3))

dev.off()

#LM
pdf(paste0(fpath,"LME_LMfrac_ZlDet_glm_gam_v2.pdf"))
par(mfrow=c(2,3))

dev.off()

#Just Zl:B GLM
pdf(paste0(fpath,"LME_All_ZlDet_glm_v2.pdf"))
par(mfrow=c(2,2))

dev.off()

#Just Zl:B GAM
pdf(paste0(fpath,"LME_All_ZlDet_gam_v2.pdf"))
par(mfrow=c(2,2))

dev.off()

#Just pTemp GLM
pdf(paste0(fpath,"LME_All_ptemp_glm_v2.pdf"))
par(mfrow=c(2,2))

dev.off()

#Just pTemp GAM
pdf(paste0(fpath,"LME_All_ptemp_gam_v2.pdf"))
par(mfrow=c(2,2))

dev.off()

#Just frac200 GLM
pdf(paste0(fpath,"LME_All_frac200_glm_v2.pdf"))
par(mfrow=c(2,2))

dev.off()

#Just frac200 GAM
pdf(paste0(fpath,"LME_All_frac200_gam_v2.pdf"))
par(mfrow=c(2,2))

dev.off()


