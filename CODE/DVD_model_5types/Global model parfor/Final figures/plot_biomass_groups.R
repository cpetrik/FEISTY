
setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor")
dataglob <- read.csv("datglob.csv",header=F)
colnames(dataglob) <- c("BioZs","BioZl","BioBs","BioBl", "BioFf","BioMf","BioPf","BioBaP","BioDf")

# load input parameters
param <- read.csv("input_parameters.csv",header=T)
outp <- cbind(param,dataglob)

##### all data loaded ###########################################################
  library(sp)
  library(rworldmap) # this pkg has waaaay better world shapefiles
  library(ggplot2)
  library(RColorBrewer)
  library(plyr)
  
  coords2 <- outp
  coords2$long <- ifelse(coords2$long > 179.5,coords2$long-360,coords2$long)
  coordinates(coords2) <- c("long", "lat")
  proj4string(coords2) <- CRS("+init=epsg:4326")
  coords2 <- spTransform(coords2, CRS("+proj=wintri"))
  test <- coordinates(coords2)
  colnames(test) <- c("long_wintri","lat_wintri")
  outp <- cbind(outp,test)

### create map
  world <- fortify(spTransform(getMap(), CRS("+proj=wintri")))
  sealand = c("#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c")
  sid = 0.15
  
  Map <- ggplot() + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=NULL),
                                    shape=15,size=sid)
  Map <- Map +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="black",fill="black")
  Map <- Map +  theme(plot.background=element_blank(),
                      panel.background=element_blank(),
                      axis.text.y   =element_blank(),
                      axis.text.x   =element_blank(),
                      axis.ticks    = element_blank(),
                      axis.title.y  =element_blank(),
                      axis.title.x  =element_blank(),
                      panel.border  = element_blank(),
                      legend.text   = element_blank(),
                      legend.title  = element_blank(),
                      legend.position = "none",
                      plot.title = element_text(hjust = 0.5))
  
  outp$BioFf <- ifelse(outp$BioFf < 10^-4, 0, outp$BioFf)
  outp$BioFf <- ifelse(outp$BioFf < 10^-2 & outp$BioFf > 0, 10^-2, outp$BioFf)
  MapBioFf <- Map + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(BioFf)),shape=15,size=sid)+
                    scale_colour_gradientn (colours= sealand,limits=c(-2,2), na.value = "white") + ggtitle("epipelagic fish")
  
  outp$BioMf <- ifelse(outp$BioMf < 10^-4, 0, outp$BioMf)
  outp$BioMf <- ifelse(outp$BioMf < 10^-2 & outp$BioMf > 0, 10^-2, outp$BioMf)
  MapBioMf <- Map + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(BioMf)),shape=15,size=sid)+
    scale_colour_gradientn (colours= sealand,limits=c(-2,2), na.value = "white") + ggtitle("mesopelagic fish")
  
  outp$BioPf <- ifelse(outp$BioPf < 10^-4, 0, outp$BioPf)
  outp$BioPf <- ifelse(outp$BioPf < 10^-2 & outp$BioPf > 0, 10^-2, outp$BioPf)
  MapBioPf <- Map + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(BioPf)),shape=15,size=sid)+
    scale_colour_gradientn (colours= sealand,limits=c(-2,2), na.value = "white") + ggtitle("large pelagic")
  
  outp$BioBaP <- ifelse(outp$BioBaP < 10^-4, 0, outp$BioBaP)
  outp$BioBaP <- ifelse(outp$BioBaP < 10^-2 & outp$BioBaP > 0, 10^-2, outp$BioBaP)
  MapBioBaP <- Map + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(BioBaP)),shape=15,size=sid)+
    scale_colour_gradientn (colours= sealand,limits=c(-2,2), na.value = "white") + ggtitle("mid-water predator")
  
  outp$BioDf <- ifelse(outp$BioDf < 10^-4, 0, outp$BioDf)
  outp$BioDf <- ifelse(outp$BioDf < 10^-2 & outp$BioDf > 0, 10^-2, outp$BioDf)
  MapBioDf <- ggplot() + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(BioDf)),shape=15,size=sid,na.rm=T)+
                         scale_colour_gradientn(colours= sealand,na.value = "white",limits=c(-2,2),
                                                labels=c("-4 to -2","-1","0","1","2"),
                         name="log10(biomass) \n (g WW m-2)")+  ggtitle("demersal fish")
  MapBioDf <- MapBioDf +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="black",fill="black")
  MapBioDf <- MapBioDf +  theme(plot.background=element_blank(),
                                panel.background=element_blank(),
                                axis.text.y   =element_blank(),
                                axis.text.x   =element_blank(),
                                axis.ticks    = element_blank(),
                                axis.title.y  =element_blank(),
                                axis.title.x  =element_blank(),
                                panel.border  = element_rect(colour = "white", size=.5,fill=NA),
                                legend.text   = element_text(size=11),
                                legend.title  = element_text(size=11),
                                plot.title = element_text(hjust = 0.5)) 

  library(cowplot)
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/Final figures/")
  jpeg(file = "Figure-sub_fish_biomass.jpeg", width=7, height=6.5,units ='in', res = 500)
  
  biomassfig <- ggdraw() + 
                draw_plot(MapBioFf,   0, .6,  .4, .3) +
                draw_plot(MapBioMf,  .4, .6,  .4, .3) +
                draw_plot(MapBioPf,   0, .3,  .4, .3) +
                draw_plot(MapBioBaP, .4, .3,  .4, .3) +
                draw_plot(MapBioDf,   0,   0, .57,  .3)

  print(biomassfig)
  dev.off()
  