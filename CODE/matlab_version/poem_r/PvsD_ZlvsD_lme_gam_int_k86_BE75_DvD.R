################################################################################

# GLMs of frac pelagic as a function of 
# Zl:D ratio, temperature, NPP, coastal area
# New parameters kt=0.0855, BE=0.075
# Use Daniel's method for GAM
# Explore interactions

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


### LINEAR MODEL ------------------------------------------------------------------
mPD <- betareg(FracPD ~ LME_ptemp + logRatZlDet + LME_Frac200 + logNPP, data=ZB)
summary(mPD) #Pseudo R-squared: 0.6421
mPDi <- betareg(FracPD ~ LME_ptemp + LME_Frac200 + logRatZlDet * logNPP, data=ZB)
summary(mPDi) #Pseudo R-squared: 0.6161
AIC(mPD,mPDi) #PDint


### ADDITIVE MODEL ------------------------------------------------------------------
## NEED TO GIVE IT BETA DISTRIBUTION 
## Daniel's methods
## NO TRANSFORMATION (never 0 or 1), PROBIT LINK (cauchit doesn't converge), 3 KNOTS

## Model selection using ZlDet
# P:D
gPD <- gam(FracPD ~ s(LME_ptemp,k=3) + s(LME_Frac200,k=3) + s(logRatZlDet,k=3) + s(logNPP,k=3), data=ZB, family=betar(link = probit))
gPDr <- gam(FracPD ~ s(logRatZlDet,k=3), data=ZB, family=betar(link = probit))
gPDt <- gam(FracPD ~ s(LME_ptemp,k=3), data=ZB, family=betar(link = probit))
gPDf <- gam(FracPD ~ s(LME_Frac200,k=3), data=ZB, family=betar(link = probit)) #edf=1 linear
gPDn <- gam(FracPD ~ s(logNPP,k=3), data=ZB, family=betar(link = probit)) 
AIC(gPD,gPDr,gPDt,gPDf,gPDn) #gPD > gPDr > gPDn
summary(gPD)

#interactions
ZB$zdnpp <- ZB$logRatZlDet * ZB$logNPP
gPDi <- gam(FracPD ~ s(LME_ptemp,k=3) + s(LME_Frac200,k=3) + s(logRatZlDet,k=3) + s(logNPP,k=3) + s(zdnpp,k=3), data=ZB, family=betar(link = probit))
AIC(gPD,gPDr,gPDt,gPDf,gPDn,gPDi) #gPDi
summary(gPDi)

# P:F
gPF <- gam(FracPF ~ s(LME_ptemp,k=3) + s(LME_Frac200,k=3) + s(logRatZlDet,k=3) + s(logNPP,k=3), data=ZB, family=betar(link = probit))
gPFr <- gam(FracPF ~ s(logRatZlDet,k=3), data=ZB, family=betar(link = probit)) #edf=1 linear
gPFt <- gam(FracPF ~ s(LME_ptemp,k=3), data=ZB, family=betar(link = probit))
gPFf <- gam(FracPF ~ s(LME_Frac200,k=3), data=ZB, family=betar(link = probit))
gPFn <- gam(FracPF ~ s(logNPP,k=3), data=ZB, family=betar(link = probit)) 
AIC(gPF,gPFr,gPFt,gPFf,gPFn) #gPF >> gPFt > gPFr > gPFn

#interactions
gPFi <- gam(FracPF ~ s(LME_ptemp,k=3) + s(LME_Frac200,k=3) + s(logRatZlDet,k=3) + s(logNPP,k=3) + s(zdnpp,k=3), data=ZB, family=betar(link = probit))
AIC(gPF,gPFr,gPFt,gPFf,gPFn,gPFi) #gPFi
summary(gPFi)

# L:M
gLM <- gam(FracLM ~ s(LME_ptemp,k=3) + s(LME_Frac200,k=3) + s(logRatZlDet,k=3) + s(logNPP,k=3), data=ZB, family=betar(link = probit))
gLMr <- gam(FracLM ~ s(logRatZlDet,k=3), data=ZB, family=betar(link = probit))
gLMt <- gam(FracLM ~ s(LME_ptemp,k=3), data=ZB, family=betar(link = probit))
gLMf <- gam(FracLM ~ s(LME_Frac200,k=3), data=ZB, family=betar(link = probit)) #edf=1 linear
gLMn <- gam(FracLM ~ s(logNPP,k=3), data=ZB, family=betar(link = probit)) 
AIC(gLM,gLMr,gLMt,gLMf,gLMn) #gLM >> gLMt >> gLMr 

#interactions
gLMi <- gam(FracLM ~ s(LME_ptemp,k=3) + s(LME_Frac200,k=3) + s(logRatZlDet,k=3) + s(logNPP,k=3) + s(zdnpp,k=3), data=ZB, family=betar(link = probit))
AIC(gLM,gLMr,gLMt,gLMf,gLMn,gLMi) #gLMi by 0.22
summary(gLMi) #interaction not sig


par(mfrow=c(3,3))
visreg(gPD)
visreg(gPDi)

par(mfrow=c(3,3))
visreg(gPF)
visreg(gPFi)

par(mfrow=c(3,3))
visreg(gLM)
visreg(gLMi)

par(mfrow=c(2,2))
visreg(gPDr)
visreg(gPFr)
visreg(gLMr)



## Model validation using ZlDet
# PD
EPD <- resid(gPDi)
FPD <- fitted(gPDi)
pdf(paste0(fpath,"LME_PDfrac_full_gam_int_resids_DvD.pdf"))
par(mfrow=c(3,3), mar=c(5,5,3,3))
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
plot(x=ZB$zdnpp,y=EPD,xlab="Zl:Det*NPP",ylab="Residuals")
abline(h=0,lty=2)
hist(EPD,xlab="",ylab="",breaks=10)
dev.off()
#no nonlinearity, just one outlier


# PF
EPF <- resid(gPFi)
FPF <- fitted(gPFi)
pdf(paste0(fpath,"LME_PFfrac_full_gam_int_resids_DvD.pdf"))
par(mfrow=c(3,3), mar=c(5,5,3,3))
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
plot(x=ZB$zdnpp,y=EPF,xlab="Zl:Det*NPP",ylab="Residuals")
abline(h=0,lty=2)
hist(EPF,xlab="",ylab="",breaks=10)
dev.off()
#one outlier


# LM
ELM <- resid(gLMi)
FLM <- fitted(gLMi)
pdf(paste0(fpath,"LME_LMfrac_full_gam_int_resids_DvD.pdf"))
par(mfrow=c(3,3), mar=c(5,5,3,3))
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
plot(x=ZB$zdnpp,y=ELM,xlab="Zl:Det*NPP",ylab="Residuals")
abline(h=0,lty=2)
hist(ELM,xlab="",ylab="",breaks=10)
dev.off()
#fitted values nonlin




### VISREG FIGURES --------------------------------------------------
#PD
pdf(paste0(fpath,"LME_PDfrac_ZlDet_gam_int_comp_DvD.pdf"))
par(mfrow=c(2,5))
visreg(gPDi,ylab="frac P vs. D",scale="response")
visreg(gPD,ylab="frac P vs. D",scale="response")
dev.off()

#PF
pdf(paste0(fpath,"LME_PFfrac_ZlDet_gam_int_comp_DvD.pdf"))
par(mfrow=c(2,5))
visreg(gPFi,ylab="frac P vs. F",scale="response")
visreg(gPF,ylab="frac P vs. F",scale="response")
dev.off()

#LM
pdf(paste0(fpath,"LME_LMfrac_ZlDet_gam_int_comp_DvD.pdf"))
par(mfrow=c(2,5))
visreg(gLMi,ylab="frac L vs. M",scale="response")
visreg(gLM,ylab="frac L vs. M",scale="response")
dev.off()



### TABLE --------------------------------------------------
library("texreg")
screenreg(list(gPDr,gPDt,gPDf,gPDn,gPD,gPDi))
screenreg(list(gPFr,gPFt,gPFf,gPFn,gPF,gPFi))
screenreg(list(gLMr,gLMt,gLMf,gLMn,gLM,gLMi))

htmlreg(list(gPDr,gPDt,gPDf,gPDn,gPD,gPDi), file = "PDgam_int_table_DvD.html")
htmlreg(list(gPFr,gPFt,gPFf,gPFn,gPF,gPFi), file = "PFgam_int_table_DvD.html")
htmlreg(list(gLMr,gLMt,gLMf,gLMn,gLM,gLMi), file = "LMgam_int_table_DvD.html")



