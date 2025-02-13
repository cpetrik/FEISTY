
rm(list=ls())

library( dplyr )

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
source("/Users/cpetrik/Dropbox/ImptDocs/R_stats/NEwR_functions_original/NumEcolR_functions/panelutils.R")
source("/Users/cpetrik/Dropbox/ImptDocs/R_stats/NEwR_functions_original/NumEcolR_functions/coldiss.R")
source("/Users/cpetrik/Dropbox/ImptDocs/R_stats/NEwR_functions_original/NumEcolR_functions/hcoplot.R")	       

setwd("/Users/cpetrik/Dropbox/Princeton/POEM_other/biol_eqs/")

# load data
smtraits <- read.csv("ICES_RAMv3.0_ALL_spp_mean_size_traits_imp_hmisc_RAMv1.0_extra.csv",sep=",",header = T,stringsAsFactors = F)

#477 stocks of 241 species
smtraits$sciname <- as.character(smtraits$sciname)
sci<-unique(smtraits$sciname)


fish <- smtraits[c("sciname","FBname","Order","EnvTemp","DemersPelag",
                   "ImpTLinf","ImpTLmax","ImpTLBoth",
                   "ImpWeight","ImpWmax","ImpWinf")]

mlo <- ddply(fish,.(Order), summarize, TL = mean( as.numeric(ImpTLBoth), na.rm = TRUE ))
mlet <- ddply(fish,.(EnvTemp), summarize, TL = mean( as.numeric(ImpTLBoth), na.rm = TRUE ))
mldp <- ddply(fish,.(DemersPelag), summarize, TL = mean( as.numeric(ImpTLBoth), na.rm = TRUE ))

mwo <- ddply(fish,.(Order), summarize, Weight = mean( as.numeric(ImpWinf), na.rm = TRUE ))
mwet <- ddply(fish,.(EnvTemp), summarize, Weight = mean( as.numeric(ImpWinf), na.rm = TRUE ))
mwdp <- ddply(fish,.(DemersPelag), summarize, Weight = mean( as.numeric(ImpWinf), na.rm = TRUE ))

mls <- ddply(fish,.(FBname), summarize, TL = mean( as.numeric(ImpTLBoth), na.rm = TRUE ))
mws <- ddply(fish,.(FBname), summarize, Weight = mean( as.numeric(ImpWinf), na.rm = TRUE ))
