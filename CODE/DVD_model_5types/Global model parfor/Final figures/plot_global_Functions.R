
library(latex2exp)

# load functions in r
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor")

  datafunc1 <- read.csv("flux_in.csv",header = F)
  datafunc1 <- datafunc1[,c(2,6,12,18,24)]
  colnames(datafunc1) <- c("FluxFf","FluxMf","FluxPf","FluxBaP","FluxDf")
  datafunc1$total <- rowSums(datafunc1,na.rm=T)
  datafunc1$potential <- datafunc1$total*0.5
  
  datafunc2 <- read.csv("foodfrac.csv",header = F)
  datafunc2$pelfrac <- datafunc2[,1]/(datafunc2[,1]+datafunc2[,2])
  
  dataglob <- read.csv("datglob.csv",header=F)
  colnames(dataglob) <- c("BioZs","BioZl","BioBs","BioBl", "BioFf","BioMf","BioPf","BioBaP","BioDf")

  param <- read.csv("input_parameters.csv",header=T)
  outp <- data.frame(param,datafunc1,pelfrac = datafunc2[,3],BioDf = dataglob$BioDf)

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

  vir <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF","#B4DE2CFF", "#FDE725FF","#ffffb2")
  sealand = c("#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c")
  sid = 0.4

  outp$potential[outp$potential < 10^-3] <- NA
  outp$potential <- ifelse(outp$potential < 10^-2, 10^-2,outp$potential)
  
  Totalflux <- ggplot() + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=log10(potential)),
                                   shape=15,size=sid)+
    scale_colour_gradientn(colours = rev(vir),limits=c(-2,2),TeX("$log_{10}$ (g WW $m^{-2}$ $y^{-1}$)"),
                                             na.value = "white") +  ggtitle("Potential fisheries production")
  Totalflux <- Totalflux +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="black",fill="black")
  Totalflux <- Totalflux +  theme(plot.background=element_blank(),
                              panel.background=element_blank(),
                              axis.text.y   =element_blank(),
                              axis.text.x   =element_blank(),
                              axis.ticks    = element_blank(),
                              axis.title.y  =element_blank(),
                              axis.title.x  =element_blank(),
                              panel.border  = element_blank(),
                              plot.title = element_text(hjust = 0.5),
                              legend.position = "bottom")

  outp$pelfrac <- ifelse(outp$BioDf < 10^-4,NA,outp$pelfrac)
  pelfrac <- ggplot() + geom_point(data=outp, aes(x=long_wintri, y=lat_wintri, colour=pelfrac),
                                     shape=15,size=sid)+
    scale_colour_gradientn (colours= sealand,limits=c(0,0.8),"Fraction of energy \n of pelagic origin",
                            na.value = "white") +  ggtitle("Benthic-pelagic coupling by demersal fish")
  pelfrac <- pelfrac +  geom_polygon(data = world, aes(x = long, y = lat, group = group),color="black",fill="black")
  pelfrac <- pelfrac +  theme(plot.background=element_blank(),
                                  panel.background=element_blank(),
                                  axis.text.y   =element_blank(),
                                  axis.text.x   =element_blank(),
                                  axis.ticks    = element_blank(),
                                  axis.title.y  =element_blank(),
                                  axis.title.x  =element_blank(),
                                  panel.border  = element_blank(),
                                  plot.title = element_text(hjust = 0.5),
                                  legend.position = "bottom")
  
  library(cowplot)
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/Final figures/")
  #pdf("Fish functions.pdf",width=12,height=4) 
  jpeg(file = "Fish functions.jpeg", width=10, height=4.5,units ='in', res = 500)
  
  
  functions <- ggdraw() + 
    draw_plot(Totalflux,    0, 0, 0.5, 1) +
    draw_plot(pelfrac,      0.5, 0, 0.5, 1) +
    draw_plot_label(c("(a)", "(b)"), c(0.02, 0.52),  c(.94, .94), size = 12, fontface = "plain")

  print(functions)
  dev.off()