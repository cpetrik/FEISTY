################################################################################

# Kiorboe & Hirst clearance data
# Adults and larvae merged
# Find types of fish (forage, pelagic, demersal)
# Fit function of weight

################################################################################

rm(list=ls())

library(lattice)
library(ggplot2)
library(ggbiplot)
library(gridExtra)

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/POEM_other/biol_eqs/")

#UAF laptop
#source(file = "/Users/Colleen/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
#setwd("/Users/Colleen/Dropbox/Princeton/POEM_other/biol_eqs/")

# load data
#clear <- read.csv("KH_clearance_merged_nohead.csv",sep=",",header = T,stringsAsFactors = F)

# smtraits <- read.csv("/Users/cpetrik/Dropbox/Princeton/RAM_Legacy/FishBase_data/RAMv1.0_mean_traits_all_orders_imputed.csv",sep=",",header = T,stringsAsFactors = F)
# obs <- c("SCIENTIFICNAME","FBname","Family","Order","Superorder","DemersPelag")
# Tr <- smtraits[obs]
# cl <- merge(clear,Tr,by.x="Species",by.y="SCIENTIFICNAME",all=T)
# cl <- cl[!is.na(cl$speciesOrig),]
# write.table(cl,"KH_clearance_merged_nohead_taxon.csv",sep=",",row.names = F)

cl <- read.csv("KH_clearance_merged_nohead_taxon.csv",sep=",",header = T,stringsAsFactors = F)
cl$Order <- as.factor(cl$Order)
cl$DemersPelag <- as.factor(cl$DemersPelag)

#----------------------------------------- PLOTS -----------------------------------------
### Specific clearance (mL / gC / day) by Order and DP
pdf("SpFmax_boxplot_Order.pdf")
Mybwplot(cl,c("mlgCdSpFmax"),"Order")
dev.off()

pdf("SpFmax_boxplot_DemersPelag.pdf")
Mybwplot(cl,c("mlgCdSpFmax"),"DemersPelag")
dev.off()

x5<-ggplot(cl, aes(x = DemersPelag, y = mlgCdSpFmax)) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1)) +
  scale_x_discrete(labels=c("benthopel","demersal","pel-neritic","reef-assoc"))
x7<-ggplot(cl, aes(x = Order, y = mlgCdSpFmax)) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1))
pdf("SpFmax_boxplot_factors.pdf")
grid.arrange(x5,x7,ncol=1)
dev.off()

x1<-ggplot(cl, aes(x = DemersPelag, y = log10(mlgCdSpFmax))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1)) +
  scale_x_discrete(labels=c("benthopel","demersal","pel-neritic","reef-assoc"))
x2<-ggplot(cl, aes(x = Order, y = log10(mlgCdSpFmax))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1))
pdf("log10SpFmax_boxplot_factors.pdf")
grid.arrange(x1,x2,ncol=1)
dev.off()

### Specific clearance (mL / gC / day) vs mass (gC)
p1 <- ggplot(cl, aes(x=MassCg, y=mlgCdSpFmax)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Mass gC") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=log10(MassCg), y=mlgCdSpFmax)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("SpFmax (mL/gC/d)")
p3 <- ggplot(cl, aes(x=log10(MassCg), y=log10(mlgCdSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("log10 SpFmax (mL/gC/d)")
pdf("SpFmax_mass.pdf")
grid.arrange(p1,p2,p3,ncol=2)
dev.off()

### Specific clearance (mL / gC / day) vs temp (C)
p1 <- ggplot(cl, aes(x=Temp, y=mlgCdSpFmax)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=Temp, y=log10(mlgCdSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("log10 SpFmax (mL/gC/d)")
pdf("SpFmax_temp.pdf")
grid.arrange(p1,p2,ncol=2)
dev.off()

#----------------------------------------- REGRESSIONS -----------------------------------------
### Mass only
m1 <- lm(mlgCdSpFmax ~ MassCg, data=cl)
summary(m1) #p-value: 0.009784

m2 <- lm(mlgCdSpFmax ~ log10(MassCg), data=cl)
summary(m2) #p-value: 0.1464

m3 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=cl)
summary(m3) #p-value: 1.978e-07

m4 <- glm(log10(mlgCdSpFmax) ~ log10(MassCg), data=cl)
summary(m4)


### See if temp effect
m5 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg) + Temp, data=cl)
summary(m5) #p-value: < 2.2e-16
anova(m5)

m6 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg) * Temp, data=cl)
summary(m6) #p-value: < 2.2e-16
anova(m6)


### See if order effect
m7 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg) + Order, data=cl)
summary(m7) #p-value: 3.723e-10
# Perc NS diff from Clupe, but Gad, Pleuro, and Salmon are
anova(m7)


### See if DP effect
m8 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg) + DemersPelag, data=cl)
summary(m8) #p-value: 3.044e-16
# Demersal NS diff from benthopelagic, but pel-neritic and reef-assoc are
anova(m8)


########################################## Q10 ###############################################
### b/c temp effect, try to adjust using assumed Q10 w/ Tref=15
cl$Q10 <- exp(0.063*(15-cl$Temp))
cl$mlgCdSpFmaxQ10 <- cl$mlgCdSpFmax * cl$Q10
cl$mlmgChSpFmaxQ10 <- cl$mlmgChSpFmax * cl$Q10

#----------------------------------------- PLOTS -----------------------------------------
### Specific clearance (mL / gC / day) by Order and DP
x1<-ggplot(cl, aes(x = DemersPelag, y = log10(mlgCdSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1)) +
  scale_x_discrete(labels=c("benthopel","demersal","pel-neritic","reef-assoc"))
x2<-ggplot(cl, aes(x = Order, y = log10(mlgCdSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1))
pdf("log10SpFmaxQ10_boxplot_factors.pdf")
grid.arrange(x1,x2,ncol=1)
dev.off()

### Specific clearance (mL / gC / day) vs mass (gC)
p1 <- ggplot(cl, aes(x=MassCg, y=mlgCdSpFmaxQ10)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Mass gC") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=log10(MassCg), y=mlgCdSpFmaxQ10)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("SpFmax (mL/gC/d)")
p3 <- ggplot(cl, aes(x=log10(MassCg), y=log10(mlgCdSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("log10 SpFmax (mL/gC/d)")
pdf("SpFmaxQ10_mass.pdf")
grid.arrange(p1,p2,p3,ncol=2)
dev.off()

### Specific clearance (mL / gC / day) vs temp (C)
p1 <- ggplot(cl, aes(x=Temp, y=mlgCdSpFmaxQ10)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=Temp, y=log10(mlgCdSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("log10 SpFmax (mL/gC/d)")
pdf("SpFmaxQ10_temp.pdf")
grid.arrange(p1,p2,ncol=2)
dev.off()

#----------------------------------------- REGRESSIONS -----------------------------------------
### Mass only
m9 <- lm(mlgCdSpFmaxQ10 ~ MassCg, data=cl)
summary(m9) #p-value: 0.004897

m10 <- lm(mlgCdSpFmaxQ10 ~ log10(MassCg), data=cl)
summary(m10) #p-value: 0.007929

m11 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=cl)
summary(m11) #p-value: 5.561e-07


### See if order effect
m12 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg) + Order, data=cl)
summary(m12) #p-value: 2.464e-09
# Gad & Salmon NS diff from Clupe, but Perc & Pleuro are
anova(m12)


### See if DP effect
m13 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg) + DemersPelag, data=cl)
summary(m13) #p-value: 5.214e-15
# Demersal NS diff from benthopelagic, but pel-neritic and reef-assoc are
anova(m13)


########################################## FORAGE ###############################################
### Simplify factors
# Clupe vs. Non-clupe
# Pel-neritic vs. non-Pel-neritic
cl$Order <- as.character(cl$Order)
cid <- which(cl$Order=="Clupeiformes")
ncid <- which(cl$Order!="Clupeiformes")
cl$ord <- cl$Order
cl$ord[ncid] <- "Other"
cl$ord <- as.factor(cl$ord)

cl$DemersPelag <- as.character(cl$DemersPelag)
pid <- which(cl$DemersPelag=="pelagic-neritic")
npid <- which(cl$DemersPelag!="pelagic-neritic")
cl$dp <- cl$DemersPelag
cl$dp[npid] <- "Other"
cl$dp <- as.factor(cl$dp)

#----------------------------------------- PLOTS -----------------------------------------
### Specific clearance (mL / gC / day) by Order and DP
x1<-ggplot(cl, aes(x = dp, y = log10(mlgCdSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1)) 
x2<-ggplot(cl, aes(x = ord, y = log10(mlgCdSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1))
pdf("log10SpFmaxQ10_boxplot_factors_simple.pdf")
grid.arrange(x1,x2,ncol=1)
dev.off()

#IN BOTH CASES, FORAGE FISHES HAVE LOWER SPECIFIC CLEARANCE THAN OTHER FISHES

#----------------------------------------- REGRESSIONS -----------------------------------------
### See if order effect
m14 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg) + ord, data=cl)
summary(m14) #p-value: 3.527e-07
# Other diff from Clupe p=0.0285
anova(m14)


### See if DP effect
m15 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg) + dp, data=cl)
summary(m15) #p-value: 2.655e-07
# Other diff from pel-neritic p=0.0206
anova(m15)

#=======================================================================================================#
#=======================================================================================================#
#=======================================================================================================#
# FIT LARVAE AND NON-LARVAE SEPARATELY
# load data
larv <- read.csv("TKClearance_rates_fish.csv",sep=",",header = T,stringsAsFactors = F)
ad <- read.csv("TKFeeding_rates_fish.csv",sep=",",header = T,stringsAsFactors = F)

#Larv same units as adults
larv$mgC <- 1e3 * larv$gC
larv$mlhClearance <- larv$mldClearance / 24

#Adjust using assumed Q10 w/ Tref=15
ad$Q10 <- exp(0.063*(15-ad$Temp))
ad$mlhFmaxQ10 <- ad$mlhFmax * ad$Q10

#Specific clearance
larv$mldClearSp <- larv$mldClearance / larv$gC
larv$mlhClearSp <- larv$mlhClearance / larv$mgC
ad$mlhFmaxSp <- ad$mlhFmaxQ10 / ad$mgC

### Specific clearance vs mass 
p1 <- ggplot(larv, aes(x=log10(mgC), y=log10(mlhClearSp))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
p2 <- ggplot(ad, aes(x=log10(mgC), y=log10(mlhFmaxSp))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
pdf("SpClear_mass_larv_ad.pdf")
grid.arrange(p1,p2,ncol=2)
dev.off()

### Clearance vs mass 
c1 <- ggplot(larv, aes(x=log10(mgC), y=log10(mlhClearance))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 Fmax (mL/h)")
c2 <- ggplot(ad, aes(x=log10(mgC), y=log10(mlhFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 Fmax (mL/h)")
pdf("Clear_mass_larv_ad.pdf")
grid.arrange(c1,c2,ncol=2)
dev.off()

### Specific clearance (mL / mgC / h) by Order and DP
x1<-ggplot(cl, aes(x = dp, y = log10(mlmgChSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1)) 
x2<-ggplot(cl, aes(x = ord, y = log10(mlmgChSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1))
pdf("log10SpFmaxQ10_boxplot_factors_simple.pdf")
grid.arrange(x1,x2,ncol=1)
dev.off()


#----------------------------------------- REGRESSIONS -----------------------------------------
### 
lc <- lm(log10(mlhClearance) ~ log10(mgC), data=larv)
summary(lc)

lsc <- lm(log10(mlhClearSp) ~ log10(mgC), data=larv)
summary(lsc)

ac <- lm(log10(mlhFmaxQ10) ~ log10(mgC), data=ad)
summary(ac)

asc <- lm(log10(mlhFmaxSp) ~ log10(mgC), data=ad)
summary(asc)


########################################## GROUPED ###############################################
Scl <- subset.data.frame(cl,MassCg<0.08)
Mcl <- subset.data.frame(cl,MassCg>=0.08 & MassCg<80)
Lcl <- subset.data.frame(cl,MassCg>80)

p1 <- ggplot(Scl, aes(x=log10(MassCmg), y=log10(mlmgChSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
p2 <- ggplot(Mcl, aes(x=log10(MassCmg), y=log10(mlmgChSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
p3 <- ggplot(cl, aes(x=log10(MassCmg), y=log10(mlmgChSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
pdf("SpClear_mass_SML.pdf")
grid.arrange(p1,p2,p3,ncol=2)
dev.off()

p1 <- ggplot(Scl, aes(x=log10(MassCmg), y=log10(mlmgChSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
p2 <- ggplot(Mcl, aes(x=log10(MassCmg), y=log10(mlmgChSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
p3 <- ggplot(cl, aes(x=log10(MassCmg), y=log10(mlmgChSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass mgC") + ylab("log10 SpFmax (mL/mgC/h)")
pdf("SpClearQ10_mass_SML.pdf")
grid.arrange(p1,p2,p3,ncol=2)
dev.off()

sc1 <- lm(log10(mlmgChSpFmax) ~ log10(MassCmg), data=Scl)
summary(sc1)

sc2 <- lm(log10(mlmgChSpFmaxQ10) ~ log10(MassCmg), data=Scl)
summary(sc2)

mc1 <- lm(log10(mlmgChSpFmax) ~ log10(MassCmg), data=Mcl)
summary(mc1)

mc2 <- lm(log10(mlmgChSpFmaxQ10) ~ log10(MassCmg), data=Mcl)
summary(mc2)

c1 <- lm(log10(mlmgChSpFmax) ~ log10(MassCmg), data=cl)
summary(c1)

c2 <- lm(log10(mlmgChSpFmaxQ10) ~ log10(MassCmg), data=cl)
summary(c2)

oc1 <- lm(log10(mlmgChSpFmax) ~ log10(MassCmg) + ord, data=cl)
summary(oc1)

oc2 <- lm(log10(mlmgChSpFmaxQ10) ~ log10(MassCmg) + ord, data=cl)
summary(oc2)

oc3 <- lm(log10(mlmgChSpFmax) ~ log10(MassCmg) * ord, data=cl)
summary(oc3)

oc4 <- lm(log10(mlmgChSpFmaxQ10) ~ log10(MassCmg) * ord, data=cl)
summary(oc4)
