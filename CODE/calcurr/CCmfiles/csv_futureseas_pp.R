################################################################################

# GLMs of recruitment var as a function of LH strat, environ, fishing
# LH classifications with k-means of 5 clusters
# Use imputed values for missing values of Fecundity & tmax
# Exclude Fecundity, PCI, and tm 
# RAM v4.4 dataset
# R/SSB
# Use cv of Fest

################################################################################

rm(list=ls())

cfile <- 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

idir <- paste0('/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/IPSL/')
gdir <- paste0('/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/GFDL/')
hdir <- paste0('/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/HAD/')

### Load data
popProdI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_pp.csv'),sep=",",header = T,stringsAsFactors = T)
popProdG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_pp.csv'),sep=",",header = T,stringsAsFactors = T)
popProdH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_pp.csv'),sep=",",header = T,stringsAsFactors = T)


all <- rbind.data.frame(popProdI,popProdG,popProdH)

dir <- paste0('/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/')
write.table(all,paste0(dir,"FEISTY_NEMURO_annual_outputs_popprod.csv"),sep=",",row.names=F)



