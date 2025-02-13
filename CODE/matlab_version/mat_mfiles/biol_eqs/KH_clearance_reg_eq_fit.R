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
cl <- read.csv("KH_clearance_merged_nohead_taxon.csv",sep=",",header = T,stringsAsFactors = F)
cl$Order <- as.factor(cl$Order)
cl$DemersPelag <- as.factor(cl$DemersPelag)

#----------------------------------------- PLOTS -----------------------------------------
### Specific clearance (mL / gC / day) by Order and DP
Mybwplot(cl,c("mlgCdSpFmax"),"Order")

Mybwplot(cl,c("mlgCdSpFmax"),"DemersPelag")

x5<-ggplot(cl, aes(x = DemersPelag, y = mlgCdSpFmax)) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1,hjust = 1)) +
  scale_x_discrete(labels=c("benthopel","demersal","pel-neritic","reef-assoc"))
x7<-ggplot(cl, aes(x = Order, y = mlgCdSpFmax)) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1,hjust = 1))
grid.arrange(x5,x7,ncol=1)

x1<-ggplot(cl, aes(x = DemersPelag, y = log10(mlgCdSpFmax))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1, 
                                                                hjust = 1)) +
  scale_x_discrete(labels=c("benthopel","demersal","pel-neritic","reef-assoc"))
x2<-ggplot(cl, aes(x = Order, y = log10(mlgCdSpFmax))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1,hjust = 1))
grid.arrange(x1,x2,ncol=1)

### Specific clearance (mL / gC / day) vs mass (gC)
p1 <- ggplot(cl, aes(x=MassCg, y=mlgCdSpFmax)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Mass gC") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=log10(MassCg), y=mlgCdSpFmax)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("SpFmax (mL/gC/d)")
p3 <- ggplot(cl, aes(x=log10(MassCg), y=log10(mlgCdSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("log10 SpFmax (mL/gC/d)")
grid.arrange(p1,p2,p3,ncol=2)

### Specific clearance (mL / gC / day) vs temp (C)
p1 <- ggplot(cl, aes(x=Temp, y=mlgCdSpFmax)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=Temp, y=log10(mlgCdSpFmax))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("log10 SpFmax (mL/gC/d)")
grid.arrange(p1,p2,ncol=2)

#----------------------------------------- REGRESSIONS -----------------------------------------
### Mass only
m1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=cl)
summary(m1) #p-value: 1.978e-07

### Simplify factors
# Clupe vs. Non-clupe
cl$Order <- as.character(cl$Order)
cid <- which(cl$Order=="Clupeiformes")
gid <- which(cl$Order=="Gadiformes")
perid <- which(cl$Order=="Perciformes")
plid <- which(cl$Order=="Pleuronectiformes")
sid <- which(cl$Order=="Salmoniformes")
cl$ord <- cl$Order
cl$ord[cid] <- "F"
cl$ord[perid] <- "P"
cl$ord[sid] <- "P"
cl$ord[gid] <- "D"
cl$ord[plid] <- "D"
cl$ord <- as.factor(cl$ord)

# Pel-neritic vs. non-Pel-neritic
cl$DemersPelag <- as.character(cl$DemersPelag)
pid <- which(cl$DemersPelag=="pelagic-neritic")
rid <- which(cl$DemersPelag=="reef-associated")
did <- which(cl$DemersPelag=="demersal")
bpid <- which(cl$DemersPelag=="benthopelagic")
cl$dp <- cl$DemersPelag
cl$dp[pid] <- "Pelag"
cl$dp[rid] <- "Pelag"
cl$dp[did] <- "Demers"
cl$dp[bpid] <- "Demers"
cl$dp <- as.factor(cl$dp)

### By Order
clup <- subset.data.frame(cl[cl$Order=="Clupeiformes",])
c1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=clup)
summary(c1) #p-value: 3.301e-07

oth <- subset.data.frame(cl[cl$Order!="Clupeiformes",])
nc1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=oth)
summary(nc1) #p-value: 0.8945

pisc <- subset.data.frame(cl[cl$ord=="P",])
p1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=pisc)
summary(p1) #p-value: 0.8397

det <- subset.data.frame(cl[cl$ord=="D",])
d1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=det)
summary(d1) #p-value: 0.7351


### By DemersPelag
pel <- subset.data.frame(cl[cl$dp=="Pelag",])
pel1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=pel)
summary(pel1) #p-value: 1.189e-08

dem <- subset.data.frame(cl[cl$dp=="Demers",])
dem1 <- lm(log10(mlgCdSpFmax) ~ log10(MassCg), data=dem)
summary(dem1) #p-value: 0.3382


########################################## Q10 ###############################################
### b/c temp effect, try to adjust using assumed Q10 w/ Tref=15
cl$Q10 <- exp(0.063*(15-cl$Temp))
cl$mlgCdSpFmaxQ10 <- cl$mlgCdSpFmax * cl$Q10
cl$mlmgChSpFmaxQ10 <- cl$mlmgChSpFmax * cl$Q10

#----------------------------------------- PLOTS -----------------------------------------
### Specific clearance (mL / gC / day) by Order and DP
x1<-ggplot(cl, aes(x = DemersPelag, y = log10(mlgCdSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1,hjust = 1)) +
  scale_x_discrete(labels=c("benthopel","demersal","pel-neritic","reef-assoc"))
x2<-ggplot(cl, aes(x = Order, y = log10(mlgCdSpFmaxQ10))) + geom_boxplot() +
  theme(text = element_text(size=12),axis.text.x = element_text(angle=45, vjust=1,hjust = 1))
grid.arrange(x1,x2,ncol=1)

### Specific clearance (mL / gC / day) vs mass (gC)
p1 <- ggplot(cl, aes(x=MassCg, y=mlgCdSpFmaxQ10)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Mass gC") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=log10(MassCg), y=mlgCdSpFmaxQ10)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("SpFmax (mL/gC/d)")
p3 <- ggplot(cl, aes(x=log10(MassCg), y=log10(mlgCdSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("log10 Mass gC") + ylab("log10 SpFmax (mL/gC/d)")
grid.arrange(p1,p2,p3,ncol=2)

### Specific clearance (mL / gC / day) vs temp (C)
p1 <- ggplot(cl, aes(x=Temp, y=mlgCdSpFmaxQ10)) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("SpFmax (mL/gC/d)")
p2 <- ggplot(cl, aes(x=Temp, y=log10(mlgCdSpFmaxQ10))) +  
  geom_point() + geom_smooth(method="lm", se=FALSE) +  xlab("Temp C") + ylab("log10 SpFmax (mL/gC/d)")
grid.arrange(p1,p2,ncol=2)

#----------------------------------------- REGRESSIONS -----------------------------------------
### Mass only
m2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=cl)
summary(m2) #p-value: 5.561e-07


### By Order
c2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=clup)
summary(c2) #p-value: 3.369e-09

nc2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=oth)
summary(nc2) #p-value: 0.5299

p2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=pisc)
summary(p2) #p-value: 0.7204

d2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=det)
summary(d2) #p-value: 0.4205


### By DemersPelag
pel2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=pel)
summary(pel2) #p-value: 0.0001407

dem2 <- lm(log10(mlgCdSpFmaxQ10) ~ log10(MassCg), data=dem)
summary(dem2) #p-value: 0.03929




