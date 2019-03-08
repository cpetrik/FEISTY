# Regression tree of parameter sensitivity test

rm(list=ls())

library(Hmisc)
#library(pastecs)
#library(fBasics)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(scatterplot3d)
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)
library(gplots)
library("corrplot")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
cfile = "Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100"
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/param_sens/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/param_sens/")

# load data
vec <- read.csv("Climatol_All_fish03_param_sens_10p_vecs_log10.csv",sep=",",header = T)

# magnitude
vec$mag <- sqrt(vec$F^2 + vec$P^2 + vec$D^2 + vec$Trop^2 + vec$Temp^2)
summary(vec$mag[2:39]) #1st Q = 0.0078

#raw data
vmat <- as.matrix(vec[2:39,c(2:6,8)])
rvec <- as.data.frame(t(vmat))
names(rvec) <- vec$Row[2:39]


## -------------------------------- correlations ----------------------------------------

# Get rid of params with small response (mag < 1st Q)
high <- which(rvec[6,]>=0.0078)
shi <- c(high,20,36)
shi <- sort(shi)
bvec <- rvec[,shi]

### Look for strong pos and neg corrs
obs1 <- names(bvec)[seq(from = 1, to = 30, by = 2)]
obs2 <- names(bvec)[seq(from = 2, to = 30, by = 2)]


pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_log10_pairplot_response_decr.pdf"))
pairs((bvec[1:5,obs1]), lower.panel = panel.cor)
dev.off()
pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_log10_pairplot_response_incr.pdf"))
pairs((bvec[1:5,obs2]), lower.panel = panel.cor)
dev.off()

tvec <- as.data.frame(t(bvec[,obs1]))

# Multiply params with corr by -1 so both changes incr F
newn <- c("ac-10","ac+10","ae-10","ae+10","am-10","am+10","alp-10","alp+10",
          "bc-10","bc+10","be-10","be+10","bm-10","bm+10","BE-10","BE+10",
          "f-10","f+10","kc-10","kc+10","ke-10","ke+10","km-10","km+10","kap-10",
          "kap+10","A-10","A+10","J-10","J+10")
names(bvec) <- newn


  
