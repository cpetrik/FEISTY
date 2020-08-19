# Fish-MIP catch & effort data for CMIP6

rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(plyr)

#PU laptop
source(file = "/Users/cpetrik/Dropbox/ImptDocs/R_stats/multiplot.r")

setwd("/Volumes/FEISTY/Fish-MIP/CMIP6/")

eff <- read.csv("Effort_PelDem_by_FGroup_AreaV4.csv",sep=",",header = T,stringsAsFactors = F)
cat <- read.csv("Catch_PelDem_by_FGroup_AreaV4.csv",sep=",",header = T,stringsAsFactors = F)


names(eff)
unique(eff$FGroup)
# "demersal<30cm"   "demersal>=90cm"  "demersal30-90cm"    
# "pelagic<30cm"    "pelagic>=90cm"   "pelagic30-90cm"

lmee <- ddply( eff, .( LMEnbr, FGroup, Year ), summarize, 
               tEffRep = sum( as.numeric(NomEffReported), na.rm = TRUE ),
               tEffIUU = sum( as.numeric(NomEffIUU), na.rm = TRUE ),
               tArea = sum( as.numeric(AreaSqKm2), na.rm = TRUE ),
               tCells = sum( as.numeric(NbrCells), na.rm = TRUE ))
lmee$tEff <- lmee$tEffIUU + lmee$tEffRep



names(cat)
unique(cat$FGroup)

lmec <- ddply( cat, .( LMEnbr, FGroup, Year ), summarize, 
               tCatRep = sum( as.numeric(Reported), na.rm = TRUE ),
               tCatIUU = sum( as.numeric(IUUs), na.rm = TRUE ),
               tArea = sum( as.numeric(AreaSqKm2), na.rm = TRUE ),
               tCells = sum( as.numeric(NbrCells), na.rm = TRUE ))
lmec$tCatch <- lmec$tCatIUU + lmec$tCatRep


write.table(lmee,"Total_effort_Fgroup_LME_year.csv",sep=",",row.names=F)
write.table(lmec,"Total_catch_Fgroup_LME_year.csv",sep=",",row.names=F)

#Combine for FEISTY functional types
# F = "pelagic<30cm"    
# P = "pelagic>=90cm"   "pelagic30-90cm"
# D = "demersal>=90cm"  "demersal30-90cm"

Flmee <- subset(lmee,FGroup=="pelagic<30cm")
Plmee <- subset(lmee,FGroup %in% c("pelagic>=90cm","pelagic30-90cm"))
Dlmee <- subset(lmee,FGroup %in% c("demersal>=90cm","demersal30-90cm"))

Flmec <- subset(lmec,FGroup=="pelagic<30cm")
Plmec <- subset(lmec,FGroup %in% c("pelagic>=90cm","pelagic30-90cm"))
Dlmec <- subset(lmec,FGroup %in% c("demersal>=90cm","demersal30-90cm"))

write.table(Flmee,"Total_effort_F_LME_year.csv",sep=",",row.names=F)
write.table(Flmec,"Total_catch_F_LME_year.csv",sep=",",row.names=F)
write.table(Plmee,"Total_effort_P_LME_year.csv",sep=",",row.names=F)
write.table(Plmec,"Total_catch_P_LME_year.csv",sep=",",row.names=F)
write.table(Dlmee,"Total_effort_D_LME_year.csv",sep=",",row.names=F)
write.table(Dlmec,"Total_catch_D_LME_year.csv",sep=",",row.names=F)

