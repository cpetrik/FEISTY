

# libraries
  library(sp)
  library(rworldmap)
  library(ggplot2)
  library(RColorBrewer)
  library(plyr)
  library(maptools)

# plot global food-web types
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor")
  dataglob <- read.csv("datglob.csv",header=F)
  colnames(dataglob) <- c("BioZs","BioZl","BioBs","BioBl", "BioFf","BioMf","BioPf","BioBaP","BioDf")
  param <- read.csv("input_parameters.csv",header=T) # load input
  outp <- cbind(param,dataglob)

  # get fraction of total biomass
  outp$tot <- outp$BioFf + outp$BioMf + outp$BioPf + outp$BioBaP + outp$BioDf
  outp$fforage <- outp$BioFf/outp$tot
  outp$fmeso <- outp$BioMf/outp$tot
  outp$fpelpred <- outp$BioPf/outp$tot
  outp$fbatpel <- outp$BioBaP/outp$tot
  outp$fdem <- outp$BioDf/outp$tot
  outp$fdempel <- outp$BioDf/(outp$BioDf+outp$BioPf)

# first shelf regions < 250 depth
  outp$group <- 0
  outp$group[outp$depthWOA <= 250 & outp$fdem >=0.05 & outp$tot  >=  10^-3] <- "shalDEM"
  outp$group[outp$group == 0 & outp$depthWOA  <= 250 & outp$fdem < 0.05 & outp$tot >= 10^-3] <- "Forage"
  outp$group[outp$group == 0 & outp$depthWOA  <= 250 & outp$tot < 10^-3] <- "None"

# now open ocean regions
  outp$group[outp$group == 0 & outp$fdem >= 0.05 & outp$tot >= 10^-3] <- "deepDEM"
  outp$group[outp$group == 0 & outp$fdem < 0.05  & outp$fpelpred >= 0.05 & outp$tot >= 10^-3] <- "PELdriv"
  outp$group[outp$group == 0 & outp$fpelpred < 0.05 & outp$fforage >= 0.05 & outp$fmeso >=0.05  & outp$tot >= 10^-3] <- "MESOFF"
  outp$group[outp$group == 0 & outp$fpelpred < 0.05 & outp$fforage < 0.05 & outp$fmeso >= 0.05  & outp$tot >= 10^-3] <- "MESO"
  outp$group[outp$group == 0 & outp$tot < 10^-3] <- "None"

# make polygons, otherwise to heavy for inkscape
  
  # assign area of interest one degrees
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
  bargrid@data$uni <- paste(tt[,1],tt[,2],sep="_")
  
  # combine with outp
  coords2 <- outp
  coords2$long2 <- ifelse(coords2$long > 179.5,coords2$long-360,coords2$long)
  coords2$uni <- paste(coords2$long2, coords2$lat ,sep="_")
  bargrid <- merge(x = bargrid, y = coords2[,c("uni" , "group")], by = "uni", all.x=TRUE)

  # select the different groups
  deepDEM <- subset(bargrid,bargrid@data$group == "deepDEM")
  Forage <- subset(bargrid,bargrid@data$group == "Forage")
  MESO <- subset(bargrid,bargrid@data$group == "MESO")
  MESOFF <- subset(bargrid,bargrid@data$group == "MESOFF")
  None <- subset(bargrid,bargrid@data$group == "None")
  PELdriv <- subset(bargrid,bargrid@data$group == "PELdriv")
  shalDEM <- subset(bargrid,bargrid@data$group == "shalDEM")
  
  deepDEM <- unionSpatialPolygons(deepDEM,deepDEM$group)
  deepDEM <- fortify(spTransform(deepDEM, CRS("+proj=wintri")))
  Forage <- unionSpatialPolygons(Forage,Forage$group)
  Forage <- fortify(spTransform(Forage, CRS("+proj=wintri")))
  MESO <- unionSpatialPolygons(MESO,MESO$group)
  MESO <- fortify(spTransform(MESO, CRS("+proj=wintri")))
  MESOFF <- unionSpatialPolygons(MESOFF,MESOFF$group)
  MESOFF <- fortify(spTransform(MESOFF, CRS("+proj=wintri")))
  None <- unionSpatialPolygons(None,None$group)
  None <- fortify(spTransform(None, CRS("+proj=wintri")))
  PELdriv <- unionSpatialPolygons(PELdriv,PELdriv$group)
  PELdriv <- fortify(spTransform(PELdriv, CRS("+proj=wintri")))
  shalDEM <- unionSpatialPolygons(shalDEM,shalDEM$group)
  shalDEM <- fortify(spTransform(shalDEM, CRS("+proj=wintri")))
  
  # get world map
  world <- fortify(spTransform(getMap(), CRS("+proj=wintri")))
  
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/Final figures/")
  pdf("Fish foodwebs world.pdf",width=7.5,height=5.5) 
  
  Mapnice <- ggplot() +  
    geom_polygon(data = Forage, aes(x = long, y = lat, group = group),color=NA,fill="#ffb380") +
    geom_polygon(data = MESO, aes(x = long, y = lat, group = group),color=NA,fill="#e0f3f8") +
    geom_polygon(data = MESOFF, aes(x = long, y = lat, group = group),color=NA,fill="#abd9e9") +
    geom_polygon(data = None, aes(x = long, y = lat, group = group),color=NA,fill="white") +
    geom_polygon(data = PELdriv, aes(x = long, y = lat, group = group),color=NA,fill="#74add1") +
    geom_polygon(data = shalDEM, aes(x = long, y = lat, group = group),color=NA,fill="#ab4300") + 
    geom_polygon(data = deepDEM2, aes(x = long, y = lat, group = group),color=NA,fill="#ff7f2a") +
    geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#808080",fill="#808080")
  
  Mapnice <- Mapnice +  theme(plot.background=element_blank(),
                              panel.background=element_blank(),
                              axis.text.y   =element_blank(),
                              axis.text.x   =element_blank(),
                              axis.title.y  =element_blank(),
                              axis.title.x  =element_blank(),
                              panel.border  = element_rect(colour = "white", size=.5,fill=NA),
                              legend.text   = element_blank(),
                              legend.title  = element_blank(),
                              legend.position = "none")
  Mapnice<- Mapnice +   guides(colour = guide_legend(override.aes = list(size=5)))
  print(Mapnice)
  dev.off()
  
  
 ########
# select locations to plot specific foodweb

#test <- subset(outp,outp$long == 190.5)
#test <- subset(test,test$lat == 20.5)
#test

#tt <- outp[18271,]
# Mapnice + geom_point(data=tt, aes(x=long_wintri, y=lat_wintri),col="red",
#                     shape=15,size=2,na.rm=T)


# north sea 739 (40 m)
# shelf break 36916 (750 m), 36272 (1100 meter), 36535 (1650 m)
# Bering Sea 18928 (50 m)
# next to Japan 12004 (3300 m) (all pels)
# peru upwelling 26410 (3800 m) (all pels)
# Chile 18271 (3800 m) (Meso-FF)
# pacific 1532 (4100 m) (Meso only)
# pacific 19713 (5200 m) (Meso only)
