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
BiomI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_biom.csv'),sep=",",header = T,stringsAsFactors = T)
RecI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_rec.csv'),sep=",",header = T,stringsAsFactors = T)
mProdI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_mprod.csv'),sep=",",header = T,stringsAsFactors = T)
tProdI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_tprod.csv'),sep=",",header = T,stringsAsFactors = T)
NcatchI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_Ncatch.csv'),sep=",",header = T,stringsAsFactors = T)
ScatchI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_Scatch.csv'),sep=",",header = T,stringsAsFactors = T)
PredI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_pred.csv'),sep=",",header = T,stringsAsFactors = T)
ecoProdI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_ecoprod.csv'),sep=",",header = T,stringsAsFactors = T)
                     

BiomG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_biom.csv'),sep=",",header = T,stringsAsFactors = T)
RecG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_rec.csv'),sep=",",header = T,stringsAsFactors = T)
mProdG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_mprod.csv'),sep=",",header = T,stringsAsFactors = T)
tProdG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_tprod.csv'),sep=",",header = T,stringsAsFactors = T)
NcatchG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_Ncatch.csv'),sep=",",header = T,stringsAsFactors = T)
ScatchG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_Scatch.csv'),sep=",",header = T,stringsAsFactors = T)
PredG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_pred.csv'),sep=",",header = T,stringsAsFactors = T)
ecoProdG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_ecoprod.csv'),sep=",",header = T,stringsAsFactors = T)


BiomH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_biom.csv'),sep=",",header = T,stringsAsFactors = T)
RecH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_rec.csv'),sep=",",header = T,stringsAsFactors = T)
mProdH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_mprod.csv'),sep=",",header = T,stringsAsFactors = T)
tProdH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_tprod.csv'),sep=",",header = T,stringsAsFactors = T)
NcatchH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_Ncatch.csv'),sep=",",header = T,stringsAsFactors = T)
ScatchH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_Scatch.csv'),sep=",",header = T,stringsAsFactors = T)
PredH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_pred.csv'),sep=",",header = T,stringsAsFactors = T)
ecoProdH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_ecoprod.csv'),sep=",",header = T,stringsAsFactors = T)


#popProdI <- read.csv(paste0(idir,'FutureSeas_IPSL_All_fishobs_pp.csv'),sep=",",header = T,stringsAsFactors = T)
#popProdG <- read.csv(paste0(gdir,'FutureSeas_GFDL_All_fishobs_pp.csv'),sep=",",header = T,stringsAsFactors = T)
#popProdH <- read.csv(paste0(hdir,'FutureSeas_HAD_All_fishobs_pp.csv'),sep=",",header = T,stringsAsFactors = T)


all <- rbind.data.frame(BiomI,RecI,mProdI,tProdI,NcatchI,ScatchI,PredI,ecoProdI,
                        BiomG,RecG,mProdG,tProdG,NcatchG,ScatchG,PredG,ecoProdG,
                        BiomH,RecH,mProdH,tProdH,NcatchH,ScatchH,PredH,ecoProdH)

dir <- paste0('/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/')
write.table(all,paste0(dir,"FEISTY_NEMURO_annual_outputs.csv"),sep=",",row.names=F)



