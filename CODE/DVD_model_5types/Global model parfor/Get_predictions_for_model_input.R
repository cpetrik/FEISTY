
rm(list=ls())

#### script to get medium and large zooplankton energy flux to fish and detritus flux to seabed
#### used as input in the fishmodel 

### input = "mg C m-2 day-1"
### output = "g C m-2 y-1"

library(R.matlab)
path <- "H:/Werk/Global drivers BP ratio/Final_analysis_manuscript/Environmental_data_Ecoregions/Charlie data" ## this is outside the git folder

#### get higher predator loss medium zooplankton 
pathname <- file.path(path, "hploss_100_mdz_1deg_ESM26_5yr_clim_191_195.mat")
datamdz <- readMat(pathname)

mdz_Month<-list(datamdz$hploss.100.mdz[1,,], datamdz$hploss.100.mdz[2,,],datamdz$hploss.100.mdz[3,,],datamdz$hploss.100.mdz[4,,],
                datamdz$hploss.100.mdz[5,,], datamdz$hploss.100.mdz[6,,],datamdz$hploss.100.mdz[7,,],datamdz$hploss.100.mdz[8,,],
                datamdz$hploss.100.mdz[9,,], datamdz$hploss.100.mdz[10,,],datamdz$hploss.100.mdz[11,,],datamdz$hploss.100.mdz[12,,])

mdz_output<-Reduce('+', mdz_Month)/12
mdz_output<-mdz_output*365/1000

#### get higher predator loss large zooplankton
pathname <- file.path(path, "hploss_100_lgz_1deg_ESM26_5yr_clim_191_195.mat")
datalgz <- readMat(pathname)

lgz_Month<-list(datalgz$hploss.100.lgz[1,,], datalgz$hploss.100.lgz[2,,],datalgz$hploss.100.lgz[3,,],datalgz$hploss.100.lgz[4,,],
                datalgz$hploss.100.lgz[5,,], datalgz$hploss.100.lgz[6,,],datalgz$hploss.100.lgz[7,,],datalgz$hploss.100.lgz[8,,],
                datalgz$hploss.100.lgz[9,,], datalgz$hploss.100.lgz[10,,],datalgz$hploss.100.lgz[11,,],datalgz$hploss.100.lgz[12,,])

lgz_output<-Reduce('+', lgz_Month)/12
lgz_output<-lgz_output*365/1000

#lgz_output[lgz_output>20] <- 20
#library(lattice)
#levelplot(lgz_output)

##### get detritus flux at seabed
pathname <- file.path(path, "fcdet_btm_1deg_ESM26_5yr_clim_191_195.mat")
datafcdet <- readMat(pathname)

fcdet_Month<-list(datafcdet$fcdet.btm[1,,], datafcdet$fcdet.btm[2,,],datafcdet$fcdet.btm[3,,],datafcdet$fcdet.btm[4,,],
                  datafcdet$fcdet.btm[5,,], datafcdet$fcdet.btm[6,,],datafcdet$fcdet.btm[7,,],datafcdet$fcdet.btm[8,,],
                  datafcdet$fcdet.btm[9,,], datafcdet$fcdet.btm[10,,],datafcdet$fcdet.btm[11,,],datafcdet$fcdet.btm[12,,])

fcdet_output<-Reduce('+', fcdet_Month)/12
fcdet_output<-fcdet_output*365/1000

##  convert to ww
lgz_output <- lgz_output*9 ## carbon to ww
mdz_output <- mdz_output*9 ## carbon to ww

fcdet_output <- fcdet_output*9 ## carbon to ww
fcdet_output <- fcdet_output*0.1 ## transfer efficiency

### get zooplankton biomass ###
pathname <- file.path(path, "mdz_100_1deg_ESM26_5yr_clim_191_195.mat")
datamdzbio <- readMat(pathname)

zooMonth<-list(datamdzbio$mdz.100[1,,], datamdzbio$mdz.100[2,,],datamdzbio$mdz.100[3,,], datamdzbio$mdz.100[4,,],
               datamdzbio$mdz.100[5,,], datamdzbio$mdz.100[6,,],datamdzbio$mdz.100[7,,], datamdzbio$mdz.100[8,,],
               datamdzbio$mdz.100[9,,], datamdzbio$mdz.100[10,,],datamdzbio$mdz.100[11,,], datamdzbio$mdz.100[12,,])
            
mdzbio_output <- Reduce('+', zooMonth)/12 # average standing stock biomass
mdzbio_output <- mdzbio_output/1000 # mg to gram
mdzbio_output <- mdzbio_output*9 # C to ww

pathname <- file.path(path, "lgz_100_1deg_ESM26_5yr_clim_191_195.mat")
datalgzbio <- readMat(pathname)

zooMonth<-list(datalgzbio$lgz.100[1,,], datalgzbio$lgz.100[2,,],datalgzbio$lgz.100[3,,], datalgzbio$lgz.100[4,,],
               datalgzbio$lgz.100[5,,], datalgzbio$lgz.100[6,,],datalgzbio$lgz.100[7,,], datalgzbio$lgz.100[8,,],
               datalgzbio$lgz.100[9,,], datalgzbio$lgz.100[10,,],datalgzbio$lgz.100[11,,], datalgzbio$lgz.100[12,,])

lgzbio_output<-Reduce('+', zooMonth)/12 #  average standing stock biomass
lgzbio_output <- lgzbio_output/1000 # mg to gram
lgzbio_output <- lgzbio_output*9 # C to ww

#### now get photic zone boundary, depth and temperature in upper 100 meters
library(raster)
setwd("H:/Werk/Global drivers BP ratio/Final_analysis_manuscript/Box model - Matcont") # outside the github folder
load("Zeu2.RData")

r_hr <- raster(nrow=4309, ncol=2149,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5)
r_hr[]<-Zeu2
r_lr <- raster(nrow=360, ncol=180,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5)
r_lr<-resample(r_hr,r_lr, method="ngb",na.rm=T)
ZeuNew<-matrix(r_lr@data@values,ncol=360,nrow=180)
ZeuNew <- cbind(ZeuNew[,181:360],ZeuNew[,1:180])

#### get depth
pathname <- file.path(path, "depth_1deg_ESM26_lat_lon.mat")
datadepth <- readMat(pathname)
pathname <- file.path(path, "depth_1deg_ESM26_lat_lon.mat")
depth <- readMat(pathname)
depth_mat <- depth$depth

### get temperature ###
pathname <- file.path(path, "temp_100_1deg_ESM26_5yr_clim_191_195.mat")
dataTemp <- readMat(pathname)

surftemp<- list(dataTemp$temp.100[1,,],dataTemp$temp.100[2,,],dataTemp$temp.100[3,,],dataTemp$temp.100[4,,],
                dataTemp$temp.100[5,,],dataTemp$temp.100[6,,],dataTemp$temp.100[7,,],dataTemp$temp.100[8,,],
                dataTemp$temp.100[9,,],dataTemp$temp.100[10,,],dataTemp$temp.100[11,,],dataTemp$temp.100[12,,])
Temp_output<-Reduce('+', surftemp)
Temp_output<-Temp_output/12

## now get datatable
lz_prod <- as.vector(lgz_output)
mz_prod <- as.vector(mdz_output)
ben_prod <- as.vector(fcdet_output)
lz_bio <- as.vector(lgzbio_output)
mz_bio <- as.vector(mdzbio_output)
photic <- as.vector(ZeuNew)
depth <- as.vector(depth_mat)
Temp_C <- as.vector(Temp_output)
lat <- rep((-89.5:89.5),360)
long <- rep((0.5:359.5),each=180)

envcond <- data.frame(long,lat,lz_prod,mz_prod,ben_prod,lz_bio,mz_bio,photic,depth,Temp_C)
envcondsub <- envcond[complete.cases(envcond[ , c(3:10)]),]

# now combine with vertical water column
pathdir <- c("H:/Werk/BP food web model shallow-deep/Ocean_temperature") # outside the github folder
setwd(pathdir)
degrees <- read.csv("woa18_A5B7_t00an01.csv",header=F,sep=",")
colnames(degrees) <- degrees[2,]
colnames(degrees)[1:3] <- c("lat","long","0")
degrees <- degrees[-c(1,2),]
degrees$lat <- as.numeric(as.character(degrees$lat))
degrees$long <- as.numeric(as.character(degrees$long))
degrees$long <- ifelse(degrees$long < 0, degrees$long + 360, degrees$long)
degrees$uni <- paste(degrees$long,degrees$lat,sep="_")
envcondsub$uni <- paste(envcondsub$long,envcondsub$lat,sep="_")
depthprof <- colnames(degrees)[3:104]
depthprof <- as.numeric(as.character(depthprof))

depthout <- c()
for (j in 1:length(envcondsub$uni)){
  nb <- subset(degrees,degrees$uni == envcondsub$uni[j])
  nb <- t(nb[,3:104])
  idx <- max(which(!is.na(nb)))
  depthout <- c(depthout, depthprof[idx])
}
envcondsub$depthWOA <- depthout 
envcondsub <- cbind(envcondsub,degrees[match(envcondsub$uni,degrees$uni),c(3:104)])
datamatlab <- envcondsub

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/")
write.csv(datamatlab, file= "input_parameters.csv")

# load output to check r.Rmax
setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Max versus realized production/")
dat_trop <- read.csv(file = "pel_prod_temp.csv",header=F)
colnames(dat_trop) <- c("FF0","FF10","FF20","FF30","PF0","PF10","PF20","PF30","ZL0","ZL10","ZL20","ZL30","ZS0","ZS10","ZS20","ZS30","bRp")
dat_trop$prodZL0   <- 1*(dat_trop$bRp - dat_trop$ZL0)
dat_trop$prodZL10  <- 1*(dat_trop$bRp - dat_trop$ZL10)
dat_trop$prodZL20  <- 1*(dat_trop$bRp - dat_trop$ZL20)
dat_trop$prodZL30  <- 1*(dat_trop$bRp - dat_trop$ZL30)

ww <- mean(dat_trop$prodZL30/ dat_trop$bRp) # 0.91 
w  <- mean(dat_trop$prodZL20/ dat_trop$bRp) # 0.89
c  <- mean(dat_trop$prodZL10/ dat_trop$bRp) # 0.85
cc <- mean(dat_trop$prodZL0/ dat_trop$bRp)  # 0.79

# use these numbers to estimate most likely RmaxL
datamatlab$RmaxL <- 0
datamatlab$RmaxL <- ifelse(datamatlab$Temp_C > 25, datamatlab$lz_prod/0.9,datamatlab$RmaxL)
datamatlab$RmaxL <- ifelse((datamatlab$Temp_C > 15 & datamatlab$Temp_C <= 25), datamatlab$lz_prod/0.89,datamatlab$RmaxL)
datamatlab$RmaxL <- ifelse((datamatlab$Temp_C > 5 & datamatlab$Temp_C <= 15), datamatlab$lz_prod/0.86,datamatlab$RmaxL)
datamatlab$RmaxL <- ifelse(datamatlab$Temp_C <=5, datamatlab$lz_prod/0.8,datamatlab$RmaxL)

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/")
write.csv(datamatlab, file= "input_parameters.csv")

datamatlab <- subset(datamatlab,!(is.na(datamatlab$depthWOA)))
datamatlab <- subset(datamatlab,!(is.na(datamatlab$X5)))

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/")
write.csv(datamatlab, file= "input_parameters.csv")

