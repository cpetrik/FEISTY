################################################################################

# GLMs of frac pelagic as a function of Z:B ratio and temperature

################################################################################

rm(list=ls())

library(visreg)

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
#cfile = "Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100"
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/")

# load data
lme <- read.csv(paste0('LME_JC15_values_',cfile,'.csv'),sep=",",header = T,stringsAsFactors = T)
lelc <- read.csv("/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/LME_Stock_PNAS_JC15_values.csv",sep=",",header = T,stringsAsFactors = T)

# Rankings of LME biomass           
corr <- cor.test(lme$JCrank, y=lme$Prank, method = 'spearman')
#S = 44870, p-value = 0.6125

corr <- cor.test(lme$JCrank, y=lme$Prank, method = 'kendall')
#z = 0.51467, p-value = 0.6068

# NPP
corr <- cor.test(lelc$JCPP, y=lelc$CNPP, method = 'pearson')
#t = 4.2269, df = 43, p-value = 0.0001212

npp<-glm(JCPP~CNPP,data=lelc)
summary(npp)
#(Intercept)   0.1287     0.1413   0.910 0.367652    
#CNPP          1.0731     0.2539   4.227 0.000121 ***

npp2<-glm(JCPP~CNPP-1,data=lelc)
summary(npp2)
#CNPP  1.29012    0.08704   14.82   <2e-16 ***

par(mfrow=c(1,2))
visreg(npp)
visreg(npp2)

t.test(lelc$JCPP,lelc$CNPP)
# t = 2.6087, df = 65.065, p-value = 0.01126
# mean of x mean of y 
# 0.6897778 0.5228925 
t.test(lelc$JCPP,lelc$CNPP, paired=TRUE)
# t = 3.4735, df = 44, p-value = 0.001165
# mean of the differences 
# 0.1668852

corr <- cor.test(lelc$JCPP, y=lelc$CNPP, method = 'spearman')
#S = 5473.4, p-value = 2.265e-06

corr <- cor.test(lelc$JCPP, y=lelc$CNPP, method = 'kendall')
#z = 4.2563, p-value = 2.078e-05



