################################################################################

# GLMs of frac pelagic as a function of 
# Zl:D ratio, temperature, NPP, coastal area
# New parameters kt=0.0855, BE=0.075

################################################################################

rm(list=ls())

library(betareg)
library(visreg)
library(lattice) #multipanel plots
library(mgcv) #gam


#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
cfile = "Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100"
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
vars <- c("LME_ptemp","LME_Frac200","logRatZlDet","logNPP")
pairs(ZB[,vars],lower.panel = panel.cor)
#lme frac200 = -0.7 with ZlDet
#npp 0.5 with ZlDet and 0.6 with lme temp
corvif(ZB[,vars])
#all<3.02

# Linear relationships?
Myxyplot(ZB,vars,"FracPD") #temp nonlin
Myxyplot(ZB,vars,"FracPF") #temp
Myxyplot(ZB,vars,"FracLM") #npp and ZlDet

### LINEAR MODEL ------------------------------------------------------------------
## FRACTION P VS. D
mPDf <- betareg(FracPD ~ LME_ptemp + logRatZlDet + LME_Frac200 + logNPP, data=ZB)
summary(mPDf) #Pseudo R-squared: 0.6421
mPDd <- betareg(FracPD ~ LME_ptemp + logRatZlDet + LME_depth + logNPP, data=ZB)
summary(mPDd) #0.6417
mPDt <- betareg(FracPD ~ LME_ptemp, data=ZB)
summary(mPDt) #0.1683
mPDr <- betareg(FracPD ~ logRatZlDet, data=ZB)
summary(mPDr) #0.5506
mPDf2 <- betareg(FracPD ~ LME_Frac200, data=ZB)
summary(mPDf2) #0.219
mPDd2 <- betareg(FracPD ~ LME_depth, data=ZB) #depth better than frac200
summary(mPDd2) #0.3519
mPDn <- betareg(FracPD ~ logNPP, data=ZB)
summary(mPDn) #0.4153

par(mfrow=c(2,2))
visreg(mPDf)
par(mfrow=c(2,2))
visreg(mPDd)
par(mfrow=c(1,1))
visreg(mPDr)
par(mfrow=c(1,2))
visreg(mPDr)
visreg(mPDn)


## FRACTION P VS. F
mPFf <- betareg(FracPF ~ LME_ptemp + logRatZlDet + LME_Frac200 + logNPP, data=ZB)
summary(mPFf)
mPFd <- betareg(FracPF ~ LME_ptemp + logRatZlDet + LME_depth + logNPP, data=ZB)
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
par(mfrow=c(1,2))
visreg(mPFr)
visreg(mPFn)


## FRACTION L VS. M
mLMf <- betareg(FracLM ~ LME_ptemp + logRatZlDet + LME_Frac200 + logNPP, data=ZB)
summary(mLMf)
mLMd <- betareg(FracLM ~ LME_ptemp + logRatZlDet + LME_depth + logNPP, data=ZB)
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
pdf(paste0(fpath,"LME_PDfrac_full_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPD,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
dev.off()
#definitely outliers

EPD <- resid(mPDr)
FPD <- fitted(mPDr)
pdf(paste0(fpath,"LME_PDfrac_ZlDet_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPD,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
dev.off()
#temp nonlinearity, outliers

# PF
EPF <- resid(mPFf)
FPF <- fitted(mPFf)
pdf(paste0(fpath,"LME_PFfrac_full_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPF,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
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
plot(x=ZB$logNPP,y=EPF,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
dev.off()
#nonlin temp

# LM
ELM <- resid(mLMf)
FLM <- fitted(mLMf)
pdf(paste0(fpath,"LME_LMfrac_full_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="LME frac area <200m",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=ELM,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
dev.off()
#def nonlinear resids with fitted

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
plot(x=ZB$logNPP,y=ELM,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
dev.off()
#resids have neg relat with temp



### ADDITIVE MODEL ------------------------------------------------------------------
## NEED TO GIVE IT BETA DISTRIBUTION (family=beta)?

## Model selection using ZlDet
gPD <- gam(FracPD ~ s(LME_ptemp,k=3) + s(logRatZlDet,k=3) + s(LME_Frac200,k=3) + s(logNPP,k=3), data=ZB, family = betar)
summary(gPD) #85.7% deviance
gPDl <- gam(FracPD ~ s(LME_ptemp,k=3) + (logRatZlDet) + s(LME_Frac200,k=3) + s(logNPP,k=3), data=ZB, family = betar)
summary(gPDl) #82.7%
gPDr <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar)
summary(gPDr) #70.9%
gPDt <- gam(FracPD ~ s(LME_ptemp,k=3), data=ZB, family = betar)
summary(gPDt) #50.2%
gPDf <- gam(FracPD ~ s(LME_Frac200,k=3), data=ZB, family = betar) #edf=1 linear
summary(gPDf) #30.4%
gPDn <- gam(FracPD ~ s(logNPP,k=3), data=ZB, family = betar) 
summary(gPDn) #54.6%
AIC(gPD,gPDr,gPDl,gPDt,gPDf,gPDn) #gPD > gPDl > gPDr > gPDn
AIC(gPD,mPDf) #gPD

gPF <- gam(FracPF ~ s(LME_ptemp,k=3) + s(logRatZlDet,k=3) + s(LME_Frac200,k=3) + s(logNPP,k=3), data=ZB, family = betar)
summary(gPF) #79.6%
gPFl <- gam(FracPF ~ s(LME_ptemp,k=3) + (logRatZlDet) + s(LME_Frac200,k=3) + s(logNPP,k=3), data=ZB, family = betar)
summary(gPFl) #75.5%
gPFr <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar) #edf=1 linear
summary(gPFr) #25.5%
gPFt <- gam(FracPF ~ s(LME_ptemp,k=3), data=ZB, family = betar)
summary(gPFt) #36.6%
gPFf <- gam(FracPF ~ s(LME_Frac200,k=3), data=ZB, family = betar)
summary(gPFf) #9.13%
gPFn <- gam(FracPF ~ s(logNPP,k=3), data=ZB, family = betar) 
summary(gPFn) #21.8%
AIC(gPF,gPFr,gPFl,gPFt,gPFf,gPFn) #gPF > gPFl >> gPFt > gPFr > gPFn
AIC(gPF,mPFf) #gPF

gLM <- gam(FracLM ~ s(LME_ptemp,k=3) + s(logRatZlDet,k=3) + s(LME_Frac200,k=3) + s(logNPP,k=3), data=ZB, family = betar)
summary(gLM) #88.1%
gLMl <- gam(FracLM ~ s(LME_ptemp,k=3) + (logRatZlDet) + s(LME_Frac200,k=3) + s(logNPP,k=3), data=ZB, family = betar)
summary(gLMl) #87.8%
gLMr <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar)
summary(gLMr) #23.9%
gLMt <- gam(FracLM ~ s(LME_ptemp,k=3), data=ZB, family = betar)
summary(gLMt) #60.5%
gLMf <- gam(FracLM ~ s(LME_Frac200,k=3), data=ZB, family = betar) #edf=1 linear
summary(gLMf) #13.2%
gLMn <- gam(FracLM ~ s(logNPP,k=3), data=ZB, family = betar) 
summary(gLMn) #20.2%
AIC(gLM,gLMr,gLMl,gLMt,gLMf,gLMn) #gLMl > gLM >> gLMt >> gLMr 
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
par(mfrow=c(2,2))
visreg(gPDl)
par(mfrow=c(2,2))
visreg(gPFl)
par(mfrow=c(2,2))
visreg(gLMl)


## Model validation using ZlDet
# PD
EPD <- resid(gPD)
FPD <- fitted(gPD)
pdf(paste0(fpath,"LME_PDfrac_full_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPD,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
dev.off()
#no nonlinearity, just one outlier

EPD <- resid(gPDr)
FPD <- fitted(gPDr)
pdf(paste0(fpath,"LME_PDfrac_ZlDet_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPD,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPD,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPD,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
dev.off()
#temp nonlinearity, outliers

# PF
EPF <- resid(gPF)
FPF <- fitted(gPF)
pdf(paste0(fpath,"LME_PFfrac_full_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPF,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
dev.off()
#one outlier

EPF <- resid(gPFr)
FPF <- fitted(gPFr)
pdf(paste0(fpath,"LME_PFfrac_ZlDet_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=EPF,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=EPF,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=EPF,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
dev.off()
#resids have nonlin relat with temp

# LM
ELM <- resid(gLM)
FLM <- fitted(gLM)
pdf(paste0(fpath,"LME_LMfrac_full_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=ELM,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
dev.off()
#fitted values nonlin

ELM <- resid(gLMr)
FLM <- fitted(gLMr)
pdf(paste0(fpath,"LME_LMfrac_ZlDet_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlDet,y=ELM,xlab="log10(Zl:Det)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_Frac200,y=ELM,xlab="Frac 200",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logNPP,y=ELM,xlab="log10(NPP)",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
dev.off()
#resids have strong relat with temp

##Test knots
gPDr4 <- gam(FracPD ~ s(logRatZlDet,k=4), data=ZB, family = betar)
gPDr5 <- gam(FracPD ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gPDr6 <- gam(FracPD ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gPDr,gPDr4,gPDr5,gPDr6) #4 > 5 > 3

gPFr4 <- gam(FracPF ~ s(logRatZlDet,k=4), data=ZB, family = betar)
gPFr5 <- gam(FracPF ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gPFr6 <- gam(FracPF ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gPFr,gPFr4,gPFr5,gPFr6) #indentical

gLMr4 <- gam(FracLM ~ s(logRatZlDet,k=4), data=ZB, family = betar)
gLMr5 <- gam(FracLM ~ s(logRatZlDet,k=5), data=ZB, family = betar)
gLMr6 <- gam(FracLM ~ s(logRatZlDet,k=6), data=ZB, family = betar)
AIC(gLMr,gLMr4,gLMr5,gLMr6) #practically identical

##Test link fun
gPDrG <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="logit"))
gPDrP <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gPDrL <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="cloglog"))
gPDrC <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="cauchit"))
AIC(gPDrG,gPDrP,gPDrL,gPDrC) #cauchit by 4

gPFrP <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gPFrL <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="cloglog"))
gPFrC <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="cauchit"))
AIC(gPFrP,gPFr,gPFrL,gPFrC) #probit>logit by 0.1

gLMrP <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="probit"))
gLMrL <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="cloglog"))
gLMrC <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family = betar(link="cauchit"))
AIC(gLMrP,gLMr,gLMrL,gLMrC) #probit by 0.005



### VISREG FIGURES --------------------------------------------------
#PD
pdf(paste0(fpath,"LME_PDfrac_ZlDet_glm_gam.pdf"))
par(mfrow=c(2,4))
visreg(mPDf,ylab="frac P vs. D")
visreg(gPD,ylab="frac P vs. D",scale="response")
dev.off()

#PF
pdf(paste0(fpath,"LME_PFfrac_ZlDet_glm_gam.pdf"))
par(mfrow=c(2,4))
visreg(mPFf,ylab="frac P vs. F")
visreg(gPF,ylab="frac P vs. F",scale="response")
dev.off()

#LM
pdf(paste0(fpath,"LME_LMfrac_ZlDet_glm_gam.pdf"))
par(mfrow=c(2,4))
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



