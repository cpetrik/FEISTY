
# estimate potential fisheries production using the flux from juv to adult

# derive size of grids
  library(rgdal)
  library(raster)
  
  gt<-(GridTopology(c(-179.5, -89.5), c(1, 1), c(360, 180))) # c(long, lat), c(cellsize long, lat), c(nb of grids long, lat)
  grt<-SpatialGrid(gt, proj4string=CRS("+init=epsg:4326"))
  spix <- as(grt, "SpatialPixels")
  spol <- as(spix, "SpatialPolygons")
  rnames<-sapply(slot(spol, "polygons"), function(x) slot(x, "ID"))
  LOCUNI<-as.data.frame(seq(1,length(spix)))
  rownames(LOCUNI)<-rnames
  bargrid<-SpatialPolygonsDataFrame(spol, LOCUNI)
  bargrid@bbox # make sure "min" is a whole number
  bargrid@data$uni <- c(1:length(bargrid))
  
  tt <- coordinates(bargrid)
  bargrid@data$long <- tt[,1]
  bargrid@data$lat <- tt[,2]
  are <- area(bargrid) # get surface area in m2
  bargrid@data$area <- are
  bargrid@data$uni <- paste(bargrid@data$long,bargrid@data$lat, sep="_")

# now load fisheries potential
setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor")

datafunc1 <- read.csv("flux_in.csv",header = F)
datafunc1 <- datafunc1[,c(2,6,12,18,24)]

colnames(datafunc1) <- c("FluxFf","FluxMf","FluxPf","FluxBaP","FluxDf")
datafunc1$potential_small <- rowSums(datafunc1[,c(1:2)],na.rm=T)
datafunc1$potential_large <- rowSums(datafunc1[,c(3:5)],na.rm=T)
datafunc1$potential_meso <- datafunc1[,c(2)]

param <- read.csv("input_parameters.csv",header=T)
outp <- data.frame(param,datafunc1)

bargrid_all <- cbind(bargrid, outp[match(bargrid@data$uni,outp$uni), c("potential_small")])
colnames(bargrid_all@data)[ncol(bargrid_all@data)] <- "potential_small"
bargrid_all <- cbind(bargrid_all, outp[match(bargrid_all@data$uni,outp$uni), c("potential_large")])
colnames(bargrid_all@data)[ncol(bargrid_all@data)] <- "potential_large"
bargrid_all <- cbind(bargrid_all, outp[match(bargrid_all@data$uni,outp$uni), c("potential_meso")])
colnames(bargrid_all@data)[ncol(bargrid_all@data)] <- "potential_meso"

bargrid_all@data$total_small <- bargrid_all@data$area * bargrid_all@data$potential_small
bargrid_all@data$total_large <- bargrid_all@data$area * bargrid_all@data$potential_large
bargrid_all@data$total_meso <- bargrid_all@data$area * bargrid_all@data$potential_meso

global_S <- sum(bargrid_all@data$total_small,na.rm =T)  # grams ww all ocean per year
global_S <- global_S / 10^6 # tonnes ww all ocean per year
global_S <- global_S / 10^6 # million tonnes ww all ocean per year
global_S

global_L <- sum(bargrid_all@data$total_large,na.rm =T)  # grams ww all ocean per year
global_L <- global_L / 10^6 # tonnes ww all ocean per year
global_L <- global_L / 10^6 # million tonnes ww all ocean per year
global_L

global_S*0.25 + global_L*0.5
global_S*0.5 + global_L*0.5

global_meso <- sum(bargrid_all@data$total_meso,na.rm =T)  # grams ww all ocean per year
global_meso <- global_meso / 10^6 # tonnes ww all ocean per year
global_meso <- global_meso / 10^6 # million tonnes ww all ocean per year
global_meso

global_meso*0.25

  # now only shelf + slope
  param <- read.csv("input_parameters.csv",header=T)
  outp_shallow <- data.frame(param,datafunc1)
  outp_shallow <- subset(outp_shallow,outp_shallow$depthWOA < 2000)
  
  bargrid <- cbind(bargrid, outp_shallow[match(bargrid@data$uni,outp_shallow$uni), c("potential_small")])
  colnames(bargrid@data)[ncol(bargrid@data)] <- "potential_small"
  bargrid <- cbind(bargrid, outp_shallow[match(bargrid@data$uni,outp_shallow$uni), c("potential_large")])
  colnames(bargrid@data)[ncol(bargrid@data)] <- "potential_large"
  bargrid <- cbind(bargrid, outp_shallow[match(bargrid@data$uni,outp_shallow$uni), c("potential_meso")])
  colnames(bargrid@data)[ncol(bargrid@data)] <- "potential_meso"
  
  bargrid@data$total_small <- bargrid@data$area * bargrid@data$potential_small
  bargrid@data$total_large <- bargrid@data$area * bargrid@data$potential_large
  bargrid@data$total_meso <- bargrid@data$area * bargrid@data$potential_meso
  
  global_S <- sum(bargrid@data$total_small,na.rm =T)  # grams ww all ocean per year
  global_S <- global_S / 10^6 # tonnes ww all ocean per year
  global_S <- global_S / 10^6 # million tonnes ww all ocean per year
  global_S
  
  global_L <- sum(bargrid@data$total_large,na.rm =T)  # grams ww all ocean per year
  global_L <- global_L / 10^6 # tonnes ww all ocean per year
  global_L <- global_L / 10^6 # million tonnes ww all ocean per year
  global_L
  
  global_S*0.25 + global_L*0.5
  global_S*0.5 + global_L*0.5
  

global_meso <- sum(bargrid@data$total_meso,na.rm =T)  # grams ww all ocean per year
global_meso <- global_meso / 10^6 # tonnes ww all ocean per year
global_meso <- global_meso / 10^6 # million tonnes ww all ocean per year
global_meso

