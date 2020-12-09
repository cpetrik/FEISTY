

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor")

# load input parameters
  param <- read.csv("input_parameters.csv",header=T)
  outp <- param

##### all data loaded ###########################################################
  library(sp)
  library(rworldmap) 
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
  valuescol <- c(0,40, 70, 80, 90 ,95, 100)/100
  sid = 0.2
  
  Depth <- ggplot() +  geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(depthWOA+1)),shape=15,size=sid)
  Depth <-  Depth +  scale_colour_gradientn (colours= rev(sealand), na.value = "white", labs(fill='log10(Depth) \n    (m)'),
                            values = valuescol)
  Depth <- Depth +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#f0f0f0",fill="#f0f0f0")
  Depth <- Depth + theme(plot.background=element_blank(),
                   panel.background=element_blank(),
                   axis.text.y   =element_blank(),
                   axis.text.x   =element_blank(),
                   axis.ticks    = element_blank(),
                   axis.title.y  =element_blank(),
                   axis.title.x  =element_blank(),
                   panel.border  = element_blank(),
                   legend.position = "bottom")
  
  outp$photic2 <- ifelse(outp$photic > 600,600,outp$photic)
  
  Phot <- ggplot() +  geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=photic2),shape=15,size=sid)
  Phot <-  Phot +  scale_colour_gradientn (colours= rev(sealand), na.value = "white", labs(fill='Zeu \n  (m)'),
                    ,breaks = c(100,350,600),labels=c("100","350",">600"))
  Phot <- Phot +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#f0f0f0",fill="#f0f0f0")
  Phot <- Phot + theme(plot.background=element_blank(),
                         panel.background=element_blank(),
                         axis.text.y   =element_blank(),
                         axis.text.x   =element_blank(),
                         axis.ticks    = element_blank(),
                         axis.title.y  =element_blank(),
                         axis.title.x  =element_blank(),
                         panel.border  = element_blank(),
                         legend.position = "bottom")
  
  Lzoop <- ggplot() +  geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(RmaxL)),shape=15,size=sid)
  Lzoop <-  Lzoop +  scale_colour_gradientn (colours= sealand, na.value = "white", 
                                             labs(fill='log10(max.large zoop. prod.)\n    (g WW m-2 y-1)'))
  Lzoop <- Lzoop +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#f0f0f0",fill="#f0f0f0")
  Lzoop <- Lzoop + theme(plot.background=element_blank(),
                       panel.background=element_blank(),
                       axis.text.y   =element_blank(),
                       axis.text.x   =element_blank(),
                       axis.ticks    = element_blank(),
                       axis.title.y  =element_blank(),
                       axis.title.x  =element_blank(),
                       panel.border  = element_blank(),
                       legend.position = "bottom")
  
  Szoop <- ggplot() +  geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(mz_prod)),shape=15,size=sid)
  Szoop <-  Szoop +  scale_colour_gradientn (colours= sealand, na.value = "white", 
                                             labs(fill='log10(max. small zoop. prod.)\n    (g WW m-2 y-1)'))
  Szoop <- Szoop +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#f0f0f0",fill="#f0f0f0")
  Szoop <- Szoop + theme(plot.background=element_blank(),
                         panel.background=element_blank(),
                         axis.text.y   =element_blank(),
                         axis.text.x   =element_blank(),
                         axis.ticks    = element_blank(),
                         axis.title.y  =element_blank(),
                         axis.title.x  =element_blank(),
                         panel.border  = element_blank(),
                         legend.position = "bottom")
  
  valuescol <- c(0,40, 70, 80, 90 ,95, 100)/100
  
  Benth <- ggplot() +  geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(ben_prod)),shape=15,size=sid)
  Benth <-  Benth +  scale_colour_gradientn (colours= sealand, na.value = "white", 
                                             labs(fill='log10(max. benthic prod.)\n     (g WW m-2 y-1)'))
  Benth <- Benth +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#f0f0f0",fill="#f0f0f0")
  Benth <- Benth + theme(plot.background=element_blank(),
                         panel.background=element_blank(),
                         axis.text.y   =element_blank(),
                         axis.text.x   =element_blank(),
                         axis.ticks    = element_blank(),
                         axis.title.y  =element_blank(),
                         axis.title.x  =element_blank(),
                         panel.border  = element_blank(),
                         legend.position = "bottom")
  
  Temp <- ggplot() +  geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=X0),shape=15,size=sid)
  Temp <-  Temp +  scale_colour_gradientn (colours= sealand, na.value = "white", 
                                             labs(fill='Surface temperature \n     (degrees celsius)'))
  Temp <- Temp +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="#f0f0f0",fill="#f0f0f0")
  Temp <- Temp + theme(plot.background=element_blank(),
                         panel.background=element_blank(),
                         axis.text.y   =element_blank(),
                         axis.text.x   =element_blank(),
                         axis.ticks    = element_blank(),
                         axis.title.y  =element_blank(),
                         axis.title.x  =element_blank(),
                         panel.border  = element_blank(),
                         legend.position = "bottom")
  
  
  library(cowplot)

  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/Final figures/")
  jpeg(file = "Figure-sub_inputpar.jpeg", width=8.5, height=10,units ='in', res = 500)
  
  inputfig <- ggdraw() + 
    draw_plot(Depth,   0, .6,  .4, .3) +
    draw_plot(Phot,  .5, .6,  .4, .3) +
    draw_plot(Lzoop,   0, .3,  .4, .3) +
    draw_plot(Szoop, .5, .3,  .4, .3) +
    draw_plot(Benth,  0, 0,  .4, .3) +
    draw_plot(Temp,  .5, 0,  .4, .3)
  
  print(inputfig)
  dev.off()
  