
### we selected four regions from EcoBase to compare with the trait-based fish model

  library(RCurl)
  library(XML)
  library(plyr)

# obtain the list of available model
  h=basicTextGatherer()
  curlPerform(url = 'http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client_3.php',writefunction=h$update)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  liste_mod<-ldply(xmlToList(data),data.frame)

# selected models two shelfs (pelagic-driven--Peru, benthic-driven North Sea)
  # and two slopes (not many are available with detailed information)
  
  # north sea shelf bounding box: (-4.680176 51,9 62)
  # ecobase model name: 457    

  # peru shelf bounding box: (-83.21095 -18.5,-66.81935 -3) ## box covers deep areas, while fish are shallow-living (only shelf included)
  # ecobase model name: 658    

  # scotland slope region bounding box: (-83.21095 -18.5,-66.81935 -3)
  # ecobase model name: 443    

  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/")
  
# north sea bounding box BOX(-4.680176 51,9 62)
  lon_min <- -4.5 --360
  lon_max <- 8.5
  lat_min <- 50.5
  lat_max <- 61.5
  
  long <- c(355.5,356.5,357.5,358.5,359.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5) # redo long as seq doesn't work
  lat  <- seq(lat_min,lat_max,1)
  coords <- merge(long,lat)
  colnames(coords) <- c("long","lat")
  coords$uni <- paste(coords$long,coords$lat,sep="_")
  coords$inout <- 100
  param <- read.csv("input_parameters.csv",header=T)
  param <- cbind(param, coords[match(param$uni,coords$uni), c("inout")])
  colnames(param)[ncol(param)] <- "inout"
  Northsea <- subset(param,param$inout == 100)
  North <- apply(Northsea[,c(5:12,14:118)], 2, FUN = median ,na.rm=T) 
  North <- matrix(North,nrow=1)
  colnames(North) <- colnames(Northsea[c(5:12,14:118)])
  North <- as.data.frame(North)
  North$depthWOA <- 100
  North[1,31:111] <- NA

######
# peru bounding box BOX(-83.21095 -18.5,-66.81935 -3)
  lon_min <- -79.5--360
  lon_max <- -66.5--360
  lat_min <- -18.5
  lat_max <- -2.5
  
  long <- seq(lon_min,lon_max,1)
  lat  <- seq(lat_min,lat_max,1)
  
  coords <- merge(long,lat)
  colnames(coords) <- c("long","lat")
  coords$uni <- paste(coords$long,coords$lat,sep="_")
  coords$inout <- 200
  param <- read.csv("input_parameters.csv",header=T)
  param <- cbind(param, coords[match(param$uni,coords$uni), c("inout")])
  colnames(param)[ncol(param)] <- "inout"
  Peru <- subset(param,param$inout == 200)
  Peru <- subset(Peru, Peru$depthWOA < 500)
  PeruN <- apply(Peru[,c(5:12,14:118)], 2, FUN = median ,na.rm=T) 
  PeruN <- matrix(PeruN,nrow=1)
  colnames(PeruN) <- colnames(Peru[c(5:12,14:118)])
  PeruN <- as.data.frame(PeruN)
  PeruN$depthWOA <- 175
  PeruN[1,34:111] <- NA

  regeco <- rbind(North,PeruN)

######
# scotland bounding box BOX(-83.21095 -18.5,-66.81935 -3)
  lon_min <- -11.5--360
  lon_max <- -3.5--360
  lat_min <- 54.5
  lat_max <- 59.5
  
  long <- seq(lon_min,lon_max,1)
  lat  <- seq(lat_min,lat_max,1)
  coords <- merge(long,lat)
  colnames(coords) <- c("long","lat")
  coords$uni <- paste(coords$long,coords$lat,sep="_")
  coords$inout <- 300
  param <- read.csv("input_parameters.csv",header=T)
  param <- cbind(param, coords[match(param$uni,coords$uni), c("inout")])
  colnames(param)[ncol(param)] <- "inout"
  Scot <- subset(param,param$inout == 300)
  ScotN <- apply(Scot[,c(5:12,14:118)], 2, FUN = median ,na.rm=T) 
  ScotN <- matrix(ScotN,nrow=1)
  colnames(ScotN) <- colnames(Scot[c(5:12,14:118)])
  ScotN <- as.data.frame(ScotN)
  ScotN$depthWOA <- 950
  ScotN[1,56:111] <- NA
  
  regeco <- rbind(regeco,ScotN)

  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots")
  write.csv(regeco,"ecopathsim_area.csv")
