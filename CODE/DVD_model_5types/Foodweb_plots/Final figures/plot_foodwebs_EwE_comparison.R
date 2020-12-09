  
## get Estimate_fluxes_per_region

############figure A food web shallow temperate ######################
 rm(list=ls())
  
  library(latex2exp)
  library(R.matlab)
  library(plotrix)

# load ecopath data
  path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/" 

  source(paste(path,"Final figures/foodwebs_EwE.R",sep="/"))

 
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/Final figures/")
  pdf("Fish foodweb EwE-Northsea.pdf",width=7.1,height=6.6) 

 # get info
  pathname <- file.path(path, "Ecopath_1.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)

  
  agg_fin <- agg_fin_northsea
  agg_fin[,3] <- NA
  agg_fin <- agg_fin[-c(2,3,5),]
  
  Zoop <- sum(biomass[1:2])
  Bent <- sum(biomass[3])
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])

  flux <- as.matrix(foodweb$mortpr)
  
  agg_fin[1,3] <- Bent
  agg_fin[2,3] <- dem
  agg_fin[3,3] <- spel+mpel+lpel+bpel
  agg_fin[4,3] <- Zoop
  
  pel <- c(5:24)
  zoop <- c(1:2)
  dem <- c(25:30)
  bent <- 3
  
  agg_fin[5,3]     <- sum(flux[pel,zoop])
  agg_fin[6,3]      <- sum(flux[pel,bent])
  agg_fin[7,3]     <- sum(flux[pel,dem])
  agg_fin[8,3]     <- sum(flux[pel,pel])
  agg_fin[9,3]     <- sum(flux[dem,zoop])
  agg_fin[10,3]     <- sum(flux[dem,bent])
  agg_fin[11,3]    <- sum(flux[dem,dem])
  agg_fin[12,3]     <- sum(flux[dem,pel])
 
  
  agg_fin[,2] <- agg_fin[,2]*1e+6 # tonnes to gram
  agg_fin[,2] <- agg_fin[,2]/1e+6 # km2 to m2
  
  
  # 
  agg_fin$bioE <- log10(agg_fin[,2]+1)
  agg_fin$biom <- log10(agg_fin[,3]+1)
  
  ecop_cex <- c(sqrt(((4*pi*1^2)*agg_fin$bioE[1:4])/ (4*pi))*5,agg_fin$bioE[5:12]*3)
  mod_cex <- c(sqrt(((4*pi*1^2)*agg_fin$biom[1:4])/ (4*pi))*5,agg_fin$biom[5:12]*3)

  colZoop <- "#2ca25f"
  colBent <- "#d95f0e"
  colPel <- "#8f0505"
  colDem <- "black"
  colMid <- "#5c88c5"
  
  
  plot(1,1,cex = ecop_cex[1], xlim = c(0,6), ylim = c(0,3),pch=16,
       col=colBent,xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  points(2,1,cex = ecop_cex[2],pch=16,col=colDem)
  points(2,2,cex = ecop_cex[3],pch=16,col=colPel)
  points(1,2,cex = ecop_cex[4],pch=16,col=colZoop)  
  
  lines(c(1,2),c(2,2),lwd=ecop_cex[5],col=colZoop)
  lines(c(1,2),c(1,2),lwd=ecop_cex[6],col=colBent)
  lines(c(1.9,1.9),c(2,1),lwd=ecop_cex[7],col=colPel)
  lines(c(1,2),c(2,1),lwd=ecop_cex[9],col=colZoop)
  lines(c(1,2),c(1,1),lwd=ecop_cex[10],col=colBent)
  lines(c(2.1,2.1),c(1,2),lwd=ecop_cex[12],col=colDem)
  draw.circle(2.1,2.1,0.2,border=colPel,lwd = ecop_cex[8])
  draw.circle(2.1,0.9,0.2,border=colDem,lwd= ecop_cex[11] )
  
  points(1,1,cex = ecop_cex[2],pch=16,col=colBent)
  points(2,1,cex = ecop_cex[2],pch=16,col=colDem)
  points(2,2,cex = ecop_cex[3],pch=16,col=colPel)
  points(1,2,cex = ecop_cex[4],pch=16,col=colZoop) 
  
  points(4,1,cex = mod_cex[2],pch=16,col=colBent)
  points(5,1,cex = mod_cex[2],pch=16,col=colDem)
  points(5,2,cex = mod_cex[3],pch=16,col=colPel)
  points(4,2,cex = mod_cex[4],pch=16,col=colZoop)  
  
  lines(c(4,5),c(2,2),lwd=mod_cex[5],col=colZoop)
  #lines(c(4,5),c(1,2),lwd=mod_cex[6],col=colBent)
  lines(c(4.9,4.9),c(2,1),lwd=mod_cex[7],col=colPel)
  lines(c(4,5),c(2,1),lwd=mod_cex[9],col=colZoop)
  lines(c(4,5),c(1,1),lwd=mod_cex[10],col=colBent)
  lines(c(5.1,5.1),c(1,2),lwd=mod_cex[12],col=colDem)
  draw.circle(5.1,2.1,0.2,border=colPel,lwd = mod_cex[8])
  draw.circle(5.1,0.9,0.2,border=colDem,lwd= mod_cex[11] )
  
  points(4,1,cex = mod_cex[2],pch=16,col=colBent)
  points(5,1,cex = mod_cex[2],pch=16,col=colDem)
  points(5,2,cex = mod_cex[3],pch=16,col=colPel)
  points(4,2,cex = mod_cex[4],pch=16,col=colZoop) 
  
  text(2.7,2.5,"North Sea",cex=2.5)
  dev.off()

  pdf("Fish foodweb EwE-Peru.pdf",width=7.1,height=6.6) 
    
  # get info peru
  pathname <- file.path(path, "Ecopath_2.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  
  agg_fin <- agg_fin_Peru
  
  agg_fin[,3] <- NA
  agg_fin <- agg_fin[-c(2,4),]
  
  Zoop <- sum(biomass[1:2])
  Bent <- sum(biomass[3])
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])
  
  flux <- as.matrix(foodweb$mortpr)
  
  agg_fin[1,3] <- Bent
  agg_fin[2,3] <- dem
  agg_fin[3,3] <- spel+mpel+lpel+bpel
  agg_fin[4,3] <- Zoop
  
  pel <- c(5:24)
  zoop <- c(1:2)
  dem <- c(25:30)
  bent <- 3
  
  agg_fin[5,3]     <- sum(flux[pel,zoop])
  agg_fin[6,3]      <- sum(flux[pel,bent])
  agg_fin[7,3]     <- sum(flux[pel,dem])
  agg_fin[8,3]     <- sum(flux[pel,pel])
  agg_fin[9,3]     <- sum(flux[dem,zoop])
  agg_fin[10,3]     <- sum(flux[dem,bent])
  agg_fin[11,3]    <- sum(flux[dem,dem])
  agg_fin[12,3]     <- sum(flux[dem,pel])
  
  # from tonnes per km2 --> to gram per m2
  agg_fin[,2] <- agg_fin[,2]*1e+6 # tonnes to gram
  agg_fin[,2] <- agg_fin[,2]/1e+6 # km2 to m2
  
  agg_fin$bioE <- log10(agg_fin[,2]+1)
  agg_fin$biom <- log10(agg_fin[,3]+1)
  
  ecop_cex <- c(sqrt(((4*pi*1^2)*agg_fin$bioE[1:4])/ (4*pi))*5,agg_fin$bioE[5:12]*3)
  mod_cex <- c(sqrt(((4*pi*1^2)*agg_fin$biom[1:4])/ (4*pi))*5,agg_fin$biom[5:12]*3)
  
  
  plot(1,1,cex = ecop_cex[1], xlim = c(0,6), ylim = c(0,3),pch=16,
       col=colBent,xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  points(2,1,cex = ecop_cex[2],pch=16,col=colDem)
  points(2,2,cex = ecop_cex[3],pch=16,col=colPel)
  points(1,2,cex = ecop_cex[4],pch=16,col=colZoop)  
  
  lines(c(1,2),c(2,2),lwd=ecop_cex[5],col=colZoop)
  #lines(c(1,2),c(1,2),lwd=ecop_cex[6],col=colBent)
  lines(c(1.9,1.9),c(2,1),lwd=ecop_cex[7],col=colPel)
  lines(c(1,2),c(2,1),lwd=ecop_cex[9],col=colZoop)
  lines(c(1,2),c(1,1),lwd=ecop_cex[10],col=colBent)
  lines(c(2.1,2.1),c(1,2),lwd=ecop_cex[12],col=colDem)
  draw.circle(2.3,2.1,0.2,border=colPel,lwd = ecop_cex[8])
  #draw.circle(2.3,0.9,0.2,border=colDem,lwd= ecop_cex[11])

  points(1,1,cex = ecop_cex[2],pch=16,col=colBent)
  points(2,1,cex = ecop_cex[2],pch=16,col=colDem)
  points(2,2,cex = ecop_cex[3],pch=16,col=colPel)
  points(1,2,cex = ecop_cex[4],pch=16,col=colZoop) 
  
  points(4,1,cex = mod_cex[2],pch=16,col=colBent)
  points(5,1,cex = mod_cex[2],pch=16,col=colDem)
  points(5,2,cex = mod_cex[3],pch=16,col=colPel)
  points(4,2,cex = mod_cex[4],pch=16,col=colZoop)  
  
  lines(c(4,5),c(2,2),lwd=mod_cex[5],col=colZoop)
  #lines(c(4,5),c(1,2),lwd=mod_cex[6],col=colBent)
  lines(c(4.9,4.9),c(2,1),lwd=mod_cex[7],col=colPel)
  lines(c(4,5),c(2,1),lwd=mod_cex[9],col=colZoop)
  lines(c(4,5),c(1,1),lwd=mod_cex[10],col=colBent)
  lines(c(5.1,5.1),c(1,2),lwd=mod_cex[12],col=colDem)
  draw.circle(5.2,2.1,0.2,border=colPel,lwd = mod_cex[8])
  draw.circle(5.2,0.9,0.2,border=colDem,lwd= mod_cex[11] )
  
  points(4,1,cex = mod_cex[2],pch=16,col=colBent)
  points(5,1,cex = mod_cex[2],pch=16,col=colDem)
  points(5,2,cex = mod_cex[3],pch=16,col=colPel)
  points(4,2,cex = mod_cex[4],pch=16,col=colZoop) 
  
  text(2.7,2.5,"Peru (shelf)",cex=2.5)
  dev.off()
  
  pdf("Fish foodweb EwE-scotland.pdf",width=7.1,height=6.6) 
  # get info scotland
  pathname <- file.path(path, "Ecopath_3.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  
  agg_fin <- agg_fin_Scot
  pel <- data.frame(Group.1 = "Pelagic",x = 0)
  agg_fin <- rbind(agg_fin[1:3,],pel,agg_fin[4:21,])
  
  agg_fin[,3] <- NA
  agg_fin <- agg_fin[-c(2,6),]
  
  Zoop <- sum(biomass[1:2])
  Bent <- sum(biomass[3])
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])
  
  flux <- as.matrix(foodweb$mortpr)
  
  agg_fin[1,3] <- Bent
  agg_fin[2,3] <- dem
  agg_fin[3,3] <- spel+lpel
  agg_fin[4,3] <-  mpel+bpel
  agg_fin[5,3] <- Zoop
    
  pel <- c(5:8,13:18)
  mid <- c(9:12,19:24)
  zoop <- c(1:2)
  dem <- c(25:30)
  bent <- 3
  
  agg_fin[6,3]     <- sum(flux[pel,zoop])
  agg_fin[7,3]      <- sum(flux[pel,bent])
  agg_fin[8,3]     <- sum(flux[pel,dem])
  agg_fin[9,3]     <- sum(flux[pel,mid])
  agg_fin[10,3]     <- sum(flux[pel,pel])
  agg_fin[11,3]     <- sum(flux[dem,zoop])
  agg_fin[12,3]     <- sum(flux[dem,bent])
  agg_fin[13,3]    <- sum(flux[dem,dem])
  agg_fin[14,3]     <- sum(flux[dem,pel])
  agg_fin[15,3]     <- sum(flux[dem,mid])
  agg_fin[16,3]     <- sum(flux[mid,zoop])
  agg_fin[17,3]     <- sum(flux[mid,bent])
  agg_fin[18,3]    <- sum(flux[mid,dem])
  agg_fin[19,3]     <- sum(flux[mid,pel])
  agg_fin[20,3]     <- sum(flux[mid,mid])
  
  # from tonnes per km2 --> to gram per m2
  agg_fin[,2] <- agg_fin[,2]*1e+6 # tonnes to gram
  agg_fin[,2] <- agg_fin[,2]/1e+6 # km2 to m2
  
  agg_fin$bioE <- log10(agg_fin[,2]+1)
  agg_fin$biom <- log10(agg_fin[,3]+1)
  
  ecop_cex <- c(sqrt(((4*pi*1^2)*agg_fin$bioE[1:5])/ (4*pi))*5,agg_fin$bioE[6:20]*3)
  mod_cex <- c(sqrt(((4*pi*1^2)*agg_fin$biom[1:5])/ (4*pi))*5,agg_fin$biom[6:20]*3)

  plot(1,1,cex = ecop_cex[1], xlim = c(0,6), ylim = c(0,4),pch=16,
       col=colBent,xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  points(2,1,cex = ecop_cex[2],pch=16,col=colDem)
  points(2,3,cex = ecop_cex[3],pch=16,col=colPel)
  points(1,3,cex = ecop_cex[4],pch=16,col=colZoop)  
  points(2,2,cex = ecop_cex[5],pch=16,col=colMid)

  # lines(c(1,2),c(3,3),lwd=ecop_cex[6],col=colZoop)
  # lines(c(1,2),c(1,3),lwd=ecop_cex[7],col=colBent)
  # lines(c(2,2),c(1,3),lwd=ecop_cex[8],col=colDem)
  # draw.circle(2.1,0.9,0.1,border=colPel,lwd= ecop_cex[9])
  # lines(c(2,2),c(2,3),lwd=ecop_cex[10],col=colMid)

  lines(c(1,2),c(3,1),lwd=ecop_cex[11],col=colZoop)
  lines(c(1,2),c(1,1),lwd=ecop_cex[12],col=colBent)
  draw.circle(2.1,0.9,0.2,border=colDem,lwd= ecop_cex[13] )
  #lines(c(2,2),c(1,3),lwd=ecop_cex[14],col=colPel)
  lines(c(1.9,1.9),c(1,2),lwd=ecop_cex[15],col=colMid)
  
  lines(c(1,2),c(3,2),lwd=ecop_cex[16],col=colZoop)
  lines(c(1,2),c(1,2),lwd=ecop_cex[17],col=colBent)
  lines(c(2.1,2.1),c(1,2),lwd=ecop_cex[18],col=colDem)
  #lines(c(2,2),c(3,2),lwd=ecop_cex[19],col=colPel) 
  draw.circle(2.3,2.0,0.2,border=colMid,lwd= ecop_cex[20])
  
  points(1,1,cex = ecop_cex[1],pch=16,col=colBent)
  points(2,1,cex = ecop_cex[2],pch=16,col=colDem)
  #points(2,3,cex = ecop_cex[3],pch=16,col=colPel)
  points(1,3,cex = ecop_cex[4],pch=16,col=colZoop) 
  points(2,2,cex = ecop_cex[5],pch=16,col=colMid) 
  
  # now add model
  points(4,1,cex = mod_cex[1],pch=16,col=colBent)
  points(5,1,cex = mod_cex[2],pch=16,col=colDem)
  points(5,3,cex = mod_cex[3],pch=16,col=colPel)
  points(4,3,cex = mod_cex[4],pch=16,col=colZoop)  
  points(5,2,cex = mod_cex[5],pch=16,col=colMid)
  
  lines(c(4,5),c(3,3),lwd=mod_cex[6],col=colZoop)
  #lines(c(4,5),c(1,3),lwd=mod_cex[7],col=colBent)
  lines(c(4.95,4.95),c(1,3),lwd=mod_cex[8],col=colDem)
  draw.circle(5.1,3.1,0.2,border=colPel,lwd= mod_cex[9])
  lines(c(4.9,4.9),c(2,3),lwd=mod_cex[10],col=colMid)
  
  lines(c(4,5),c(3,1),lwd=mod_cex[11],col=colZoop)
  lines(c(4,5),c(1,1),lwd=mod_cex[12],col=colBent)
  draw.circle(5.1,0.9,0.2,border=colDem,lwd= mod_cex[13] )
  lines(c(5.05,5.05),c(1,3),lwd=mod_cex[14],col=colPel)
  lines(c(5,5),c(1,2),lwd=mod_cex[15],col=colMid)
  
  lines(c(4,5),c(3,2),lwd=mod_cex[16],col=colZoop)
  #lines(c(4,5),c(1,2),lwd=mod_cex[17],col=colBent)
  lines(c(5.1,5.1),c(1,2),lwd=mod_cex[18],col=colDem)
  lines(c(5,5),c(3,2),lwd=mod_cex[19],col=colPel) 
  draw.circle(5.3,2,0.2,border=colMid,lwd= mod_cex[20])
  
  points(4,1,cex = mod_cex[1],pch=16,col=colBent)
  points(5,1,cex = mod_cex[2],pch=16,col=colDem)
  points(5,3,cex = mod_cex[3],pch=16,col=colPel)
  points(4,3,cex = mod_cex[4],pch=16,col=colZoop) 
  points(5,2,cex = mod_cex[5],pch=16,col=colMid) 
  
  text(2.7,3.5,"Scotland (slope)",cex=2.5)
  dev.off()
  