################################################################################

# GLMs of frac pelagic as a function of Z:B ratio and temperature

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
ZB <- read.csv("LME_ZBratios_clim_fished_All_fish03.csv",sep=",",header = T,stringsAsFactors = T)
ZB$logRatZDet <- log10(ZB$RatZDet)
ZB$logRatZB <- log10(ZB$RatZB)
ZB$logRatZlDet <- log10(ZB$RatZlDet)
ZB$logRatZlB <- log10(ZB$RatZlB)

### DATA EXPLORATION ------------------------------------------------------------------
# Outliers
Mydotplot(as.matrix(ZB))
#RatZDet has 2: ~0.75, 2.25
#RatZB has 1 >> than the rest
#RatZlDet has 1 ~5

# Collinearity
pairs(ZB,lower.panel = panel.cor)

# Linear relationships?
vars <- c("LME_ptemp","logRatZDet","logRatZB","logRatZlDet","logRatZlB")
Myxyplot(ZB,vars,"FracPD") #temp, ZDet
Myxyplot(ZB,vars,"FracPF") #temp, ZDet
Myxyplot(ZB,vars,"FracLM") #all ratios?

### LINEAR MODEL ------------------------------------------------------------------
## FRACTION P VS. D
mPD1 <- betareg(FracPD ~ LME_ptemp + logRatZDet, data=ZB)
summary(mPD1)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     0.80255    0.48285   1.662   0.0965 .  
# LME_ptemp       0.01569    0.01481   1.060   0.2893    
# log10(RatZDet)  1.92107    0.47442   4.049  5.14e-05 ***
#Pseudo R-squared: 0.08525

mPD1i <- betareg(FracPD ~ LME_ptemp * logRatZDet, data=ZB)
summary(mPD1i) #int sig

mPD2 <- betareg(FracPD ~ LME_ptemp + logRatZB, data=ZB)
summary(mPD2)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  25.81388    3.85453   6.697  2.13e-11 ***
# LME_ptemp     0.04299    0.01403   3.063  0.00219 ** 
# log10(RatZB)  4.06902    0.58504   6.955  3.52e-12 ***
#Pseudo R-squared: 0.2831

mPD2i <- betareg(FracPD ~ LME_ptemp * logRatZB, data=ZB)
summary(mPD2i) #int sig

mPD3 <- betareg(FracPD ~ LME_ptemp + logRatZlDet, data=ZB)
summary(mPD3)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     -0.34337    0.25397  -1.352   0.1764    
# LME_ptemp       -0.03090    0.01485  -2.081   0.0375 *  
# log10(RatZlDet)  3.12740    0.49095   6.370   1.89e-10 ***
#Pseudo R-squared: 0.2397

mPD3i <- betareg(FracPD ~ LME_ptemp * logRatZlDet, data=ZB)
summary(mPD3i)

mPD4 <- betareg(FracPD ~ LME_ptemp + logRatZlB, data=ZB)
summary(mPD4)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   28.83646    2.59612  11.108  < 2e-16 ***
# LME_ptemp     -0.04621    0.01250  -3.697  0.000218 ***
# log10(RatZlB)  5.07584    0.44297  11.459  < 2e-16 ***
#Pseudo R-squared: 0.3799

mPD4i <- betareg(FracPD ~ LME_ptemp * logRatZlB, data=ZB)
summary(mPD4i) 

mPDt <- betareg(FracPD ~ LME_ptemp, data=ZB)
summary(mPDt) 

par(mfrow=c(1,2))
visreg(mPD1)
par(mfrow=c(1,2))
visreg(mPD2)
par(mfrow=c(1,2))
visreg(mPD3)
par(mfrow=c(1,2))
visreg(mPD4)


## FRACTION P VS. F
mPF1 <- betareg(FracPF ~ LME_ptemp + logRatZDet, data=ZB)
summary(mPF1)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)     0.71634    0.49486   1.448   0.1477  
# LME_ptemp      -0.02744    0.01542  -1.779   0.0752 .
# logRatZDet)  1.06359    0.47264   2.250   0.0244 *
# Pseudo R-squared: 0.08231

mPF1i <- betareg(FracPF ~ LME_ptemp * logRatZDet, data=ZB)
summary(mPF1i)

mPF2 <- betareg(FracPF ~ LME_ptemp + logRatZB, data=ZB)
summary(mPF2)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  18.895553   3.937484   4.799 1.60e-06 ***
# LME_ptemp    -0.009806   0.014763  -0.664 0.507    
# logRatZB)  2.895888   0.594508   4.871 1.11e-06 ***
# Pseudo R-squared: 0.2364

mPF2i <- betareg(FracPF ~ LME_ptemp * logRatZB, data=ZB)
summary(mPF2i)

mPF3 <- betareg(FracPF ~ LME_ptemp + logRatZlDet, data=ZB)
summary(mPF3)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      0.30936    0.26321   1.175 0.24    
# LME_ptemp       -0.06586    0.01569  -4.197 2.70e-05 ***
# logRatZlDet)  2.45503    0.49652   4.944 7.63e-07 ***
# Pseudo R-squared: 0.2151

mPF3i <- betareg(FracPF ~ LME_ptemp * logRatZlDet, data=ZB)
summary(mPF3i)

mPF4 <- betareg(FracPF ~ LME_ptemp + logRatZlB, data=ZB)
summary(mPF4)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   24.78232    2.76151   8.974   < 2e-16 ***
# LME_ptemp     -0.08929    0.01417  -6.302   2.94e-10 ***
# logRatZlB)  4.22590    0.46542   9.080   < 2e-16 ***
# Pseudo R-squared: 0.3657

mPF4i <- betareg(FracPF ~ LME_ptemp * logRatZlB, data=ZB)
summary(mPF4i) #int sig

mPFt <- betareg(FracPF ~ LME_ptemp, data=ZB)
summary(mPFt)

AIC(mPF4,mPF4i) #interaction is better

par(mfrow=c(1,2))
visreg(mPF1)
par(mfrow=c(1,2))
visreg(mPF2)
par(mfrow=c(1,2))
visreg(mPF3)
par(mfrow=c(1,2))
visreg(mPF4)
par(mfrow=c(1,2))
visreg(mPF4i)


## FRACTION L VS. M
mLM1 <- betareg(FracLM ~ LME_ptemp + logRatZDet, data=ZB)
summary(mLM1)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.167109   0.169249  -0.987   0.323    
# LME_ptemp      -0.073304   0.005551 -13.204   < 2e-16 ***
# logRatZDet) -1.169251   0.164848  -7.093   1.31e-12 ***
# Pseudo R-squared: 0.7579

mLM1i <- betareg(FracLM ~ LME_ptemp * logRatZDet, data=ZB)
summary(mLM1i)

mLM2 <- betareg(FracLM ~ LME_ptemp + logRatZB, data=ZB)
summary(mLM2)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.42145    1.67587  -2.638  0.00833 ** 
# LME_ptemp    -0.06855    0.00685 -10.007  < 2e-16 ***
# logRatZB) -0.79604    0.25319  -3.144  0.00167 ** 
# Pseudo R-squared: 0.6107

mLM2i <- betareg(FracLM ~ LME_ptemp * logRatZB, data=ZB)
summary(mLM2i) #int sig

mLM3 <- betareg(FracLM ~ LME_ptemp + logRatZlDet, data=ZB)
summary(mLM3)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      0.75803    0.12226   6.200   5.64e-10 ***
# LME_ptemp       -0.05751    0.00717  -8.021   1.05e-15 ***
# logRatZlDet) -0.47362    0.21975  -2.155   0.0311 *  
# Pseudo R-squared: 0.5824

mLM3i <- betareg(FracLM ~ LME_ptemp * logRatZlDet, data=ZB)
summary(mLM3i) #int sig

mLM4 <- betareg(FracLM ~ LME_ptemp + logRatZlB, data=ZB)
summary(mLM4)
# Coefficients (mean model with logit link):
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    2.693703   1.301767   2.069   0.0385 *  
# LME_ptemp     -0.067195   0.007458  -9.010   <2e-16 ***
# logRatZlB)  0.309291   0.218174   1.418   0.1563    
# Pseudo R-squared: 0.5644

mLM4i <- betareg(FracLM ~ LME_ptemp * logRatZlB, data=ZB)
summary(mLM4i)

mLMt <- betareg(FracLM ~ LME_ptemp, data=ZB)
summary(mLMt)

AIC(mLM1,mLMt)
AIC(mLM2,mLMt)
AIC(mLM3,mLMt)
AIC(mLM4,mLM4i,mLMt)

par(mfrow=c(1,2))
visreg(mLM1)
par(mfrow=c(1,2))
visreg(mLM2)
par(mfrow=c(1,2))
visreg(mLM3)
par(mfrow=c(1,2))
visreg(mLM4)
par(mfrow=c(1,1))
visreg(mLMt)


## Model selection using ZlB
mPD4 <- betareg(FracPD ~ LME_ptemp + logRatZlB, data=ZB)
summary(mPD4)
mPD4i <- betareg(FracPD ~ LME_ptemp * logRatZlB, data=ZB)
summary(mPD4i)
mPD4r <- betareg(FracPD ~ logRatZlB, data=ZB)
summary(mPD4r)
AIC(mPD4,mPD4i,mPD4r) #mPD4

mPF4 <- betareg(FracPF ~ LME_ptemp + logRatZlB, data=ZB)
summary(mPF4)
mPF4i <- betareg(FracPF ~ LME_ptemp * logRatZlB, data=ZB)
summary(mPF4i)
mPF4r <- betareg(FracPF ~ logRatZlB, data=ZB)
summary(mPF4r)
AIC(mPF4,mPF4i,mPF4r) #mPF4i

mLM4 <- betareg(FracLM ~ LME_ptemp + logRatZlB, data=ZB)
summary(mLM4)
mLM4i <- betareg(FracLM ~ LME_ptemp * logRatZlB, data=ZB)
summary(mLM4i)
mLM4r <- betareg(FracLM ~ logRatZlB, data=ZB)
summary(mLM4r)
AIC(mLM4,mLM4i,mLM4r,mLMt) #mLMt


par(mfrow=c(2,2))
visreg(mPD4r)
visreg(mPF4r)
visreg(mLM4r)


## Model validation using ZlB
# PD
EPD <- resid(mPD4)
FPD <- fitted(mPD4)
pdf(paste0(fpath,"LME_PDfrac_ZlB_ptemp_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#looks like some nonlinearity in the residuals and definitely outliers

EPD <- resid(mPD4r)
FPD <- fitted(mPD4r)
pdf(paste0(fpath,"LME_PDfrac_ZlB_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#no nonlinearity, just outliers

# PF
EPF <- resid(mPF4)
FPF <- fitted(mPF4)
pdf(paste0(fpath,"LME_PFfrac_ZlB_ptemp_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
# Doesn't look as bad as PD

EPF <- resid(mPF4r)
FPF <- fitted(mPF4r)
pdf(paste0(fpath,"LME_PFfrac_ZlB_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have neg relat with temp

# LM
ELM <- resid(mLMt)
FLM <- fitted(mLMt)
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
#nonlinear residuals

ELM <- resid(mLM4)
FLM <- fitted(mLM4)
pdf(paste0(fpath,"LME_LMfrac_ZlB_ptemp_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()

ELM <- resid(mLM4r)
FLM <- fitted(mLM4r)
pdf(paste0(fpath,"LME_LMfrac_ZlB_glm_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have neg relat with temp



### ADDITIVE MODEL ------------------------------------------------------------------
## NEED TO GIVE IT BETA DISTRIBUTION (family=beta)?

## Model selection using ZlB
gPD4 <- gam(FracPD ~ s(LME_ptemp,k=4) + s(logRatZlB,k=4), data=ZB, family = betar)
summary(gPD4)
gPD4r <- gam(FracPD ~ s(logRatZlB,k=4), data=ZB, family = betar)
summary(gPD4r) #82.4% dev
gPD4l <- gam(FracPD ~ s(LME_ptemp,k=4) + (logRatZlB), data=ZB, family = betar)
summary(gPD4l)
gPDt <- gam(FracPD ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gPDt)
AIC(gPD4,gPD4r,gPD4l,gPDt) #gPD4
AIC(gPD4,mPD4) #mPD4

gPF4 <- gam(FracPF ~ s(LME_ptemp,k=4) + s(logRatZlB,k=4), data=ZB, family = betar)
summary(gPF4)
gPF4r <- gam(FracPF ~ s(logRatZlB,k=4), data=ZB, family = betar)
summary(gPF4r) #53.7% dev
gPF4l <- gam(FracPF ~ s(LME_ptemp,k=4) + (logRatZlB), data=ZB, family = betar)
summary(gPF4l)
gPFt <- gam(FracPF ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gPFt)
AIC(gPF4,gPF4r,gPF4l,gPFt) #gPF4
AIC(gPF4,mPF4) #mPF4

gLM4 <- gam(FracLM ~ s(LME_ptemp,k=4) + s(logRatZlB,k=4), data=ZB, family = betar)
summary(gLM4)
gLM4r <- gam(FracLM ~ s(logRatZlB,k=4), data=ZB, family = betar)
summary(gLM4r) #17.7% dev
gLM4l <- gam(FracLM ~ s(LME_ptemp,k=4) + (logRatZlB), data=ZB, family = betar)
summary(gLM4l)
gLMt <- gam(FracLM ~ s(LME_ptemp,k=4), data=ZB, family = betar)
summary(gLMt)
AIC(gLM4,gLM4r,gLM4l,gLMt) #gLM4

par(mfrow=c(1,2))
visreg(gPD4)
par(mfrow=c(1,2))
visreg(gPF4)
par(mfrow=c(1,2))
visreg(gLM4)
#All have a kind of quadratic, concave-down relat with temp

par(mfrow=c(2,2))
visreg(gPD4r)
visreg(gPF4r)
visreg(gLM4r)
#I think Zl:B rat could be linear

#Zl:B rat linear
par(mfrow=c(1,2))
visreg(gPD4l)
visreg(gPF4l)
visreg(gLM4l)
#All have a kind of quadratic, concave-down relat with temp


## Model validation using ZlB
# PD
EPD <- resid(gPD4)
FPD <- fitted(gPD4)
pdf(paste0(fpath,"LME_PDfrac_ZlB_ptemp_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#looks okay

EPD <- resid(gPD4r)
FPD <- fitted(gPD4r)
pdf(paste0(fpath,"LME_PDfrac_ZlB_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPD,y=EPD,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPD,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPD,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
plot(sort(EPD),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#no nonlinearity, just outliers

# PF
EPF <- resid(gPF4)
FPF <- fitted(gPF4)
pdf(paste0(fpath,"LME_PFfrac_ZlB_ptemp_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#OK

EPF <- resid(gPF4r)
FPF <- fitted(gPF4r)
pdf(paste0(fpath,"LME_PFfrac_ZlB_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FPF,y=EPF,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=EPF,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=EPF,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
plot(sort(EPF),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have strong neg relat with temp

# LM
ELM <- resid(gLM4)
FLM <- fitted(gLM4)
pdf(paste0(fpath,"LME_LMfrac_ZlB_ptemp_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#OK

ELM <- resid(gLM4r)
FLM <- fitted(gLM4r)
pdf(paste0(fpath,"LME_LMfrac_ZlB_gam_resids.pdf"))
par(mfrow=c(2,3), mar=c(5,5,3,3))
plot(x=FLM,y=ELM,xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$logRatZlB,y=ELM,xlab="log10(Zl:B)",ylab="Residuals")
abline(h=0,lty=2)
plot(x=ZB$LME_ptemp,y=ELM,xlab="LME temp",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
plot(sort(ELM),type="h",xlab="Sorted residuals",ylab="Residuals")
dev.off()
#resids have strong neg relat with temp



### FIGURES --------------------------------------------------
#PD
pdf(paste0(fpath,"LME_PDfrac_ZlB_glm_gam.pdf"))
par(mfrow=c(2,2))
visreg(mPD4,ylab="frac P vs. D")
visreg(gPD4,ylab="frac P vs. D")
dev.off()

#PF
pdf(paste0(fpath,"LME_PFfrac_ZlB_glm_gam.pdf"))
par(mfrow=c(2,2))
visreg(mPF4,ylab="frac P vs. F")
visreg(gPF4,ylab="frac P vs. F")
dev.off()

#LM
pdf(paste0(fpath,"LME_LMfrac_ZlB_glm_gam.pdf"))
par(mfrow=c(2,2))
visreg(mLM4,ylab="frac L vs. M")
visreg(gLM4,ylab="frac L vs. M")
dev.off()

#Just Zl:B GLM
pdf(paste0(fpath,"LME_All_ZlB_glm.pdf"))
par(mfrow=c(2,2))
visreg(mPD4r,ylab="frac P vs. D")
visreg(mPF4r,ylab="frac P vs. F")
visreg(mLM4r,ylab="frac L vs. M")
dev.off()

#Just Zl:B GAM
pdf(paste0(fpath,"LME_All_ZlB_gam.pdf"))
par(mfrow=c(2,2))
visreg(gPD4r,ylab="frac P vs. D")
visreg(gPF4r,ylab="frac P vs. F")
visreg(gLM4r,ylab="frac L vs. M")
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
visreg(gPDt,ylab="frac P vs. D")
visreg(gPFt,ylab="frac P vs. F")
visreg(gLMt,ylab="frac L vs. M")
dev.off()

