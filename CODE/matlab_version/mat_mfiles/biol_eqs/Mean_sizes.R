
rm(list=ls())

library( plyr )

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
source("/Users/cpetrik/Dropbox/ImptDocs/R_stats/NEwR_functions_original/NumEcolR_functions/panelutils.R")
source("/Users/cpetrik/Dropbox/ImptDocs/R_stats/NEwR_functions_original/NumEcolR_functions/coldiss.R")
source("/Users/cpetrik/Dropbox/ImptDocs/R_stats/NEwR_functions_original/NumEcolR_functions/hcoplot.R")	       

setwd("/Users/cpetrik/Dropbox/Princeton/RAM_Legacy/")

# load data
smtraits <- read.csv("FishBase_data/v1_RAM/RAMv1.0_mean_traits_all_stocks_orders_PCscore_imputed.csv",sep=",",header = T,stringsAsFactors = F)
# My matches
matchstks <- read.csv("Traits_Recruit/RAM_v1/RAMv1.0_allRAMstocks_matched_FB_stock_means.csv",sep=",", header = T, stringsAsFactors = F)

#502 stocks of 310 species
smtraits$sciname <- as.character(smtraits$sciname)
sci<-unique(smtraits$sciname)

length <- c("TL","TLcom","TLMatMin","TLMatMin2","Lm","TLmax","TLinf",
            "TLFecMin","TLFecMax","TLMin","TLMax")
weight <- c("Weight","Wmax","Winf","WeightMin","WeightMax")

summary(smtraits[,weight]) #Weight
summary(smtraits[,length]) #TL

fish <- smtraits[c("sciname","FBname","Order","EnvTemp","DemersPelag","Weight","TL")]

mlo <- ddply(fish,.(Order), summarize, TL = mean( as.numeric(TL), na.rm = TRUE ))
mlet <- ddply(fish,.(EnvTemp), summarize, TL = mean( as.numeric(TL), na.rm = TRUE ))
mldp <- ddply(fish,.(DemersPelag), summarize, TL = mean( as.numeric(TL), na.rm = TRUE ))

mwo <- ddply(fish,.(Order), summarize, Weight = mean( as.numeric(Weight), na.rm = TRUE ))
mwet <- ddply(fish,.(EnvTemp), summarize, Weight = mean( as.numeric(Weight), na.rm = TRUE ))
mwdp <- ddply(fish,.(DemersPelag), summarize, Weight = mean( as.numeric(Weight), na.rm = TRUE ))

mls <- ddply(fish,.(FBname), summarize, TL = mean( as.numeric(TL), na.rm = TRUE ))
mws <- ddply(fish,.(FBname), summarize, Weight = mean( as.numeric(Weight), na.rm = TRUE ))
