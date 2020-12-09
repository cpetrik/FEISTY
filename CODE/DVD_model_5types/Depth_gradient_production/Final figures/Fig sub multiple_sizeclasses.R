rm(list=ls())
library(latex2exp)

#### script to get matlab biomass output 

library(R.matlab)
path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Depth_gradient_production"

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Depth_gradient_production/Final figures")
pdf("Figure_S21.pdf",width=5,height=6.5)
op <- par(mfrow = c(4,3),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

linesize  <- 1.5

###################################################################################
#### get biomass
pathname <- file.path(path, "depth_threeclasses.mat")
databio <- readMat(pathname)

forage  <- databio$BioFf
meso    <- databio$BioMf
pelpred <- databio$BioPf
batpel  <- databio$BioBaP
dem     <- databio$BioDf

depth <- as.numeric(databio$bdepth)
depth <- log10(depth)
pelprod <- as.numeric(databio$bRp)
benprod <- as.numeric(databio$bRb)
ph <- databio$param[,,1]$photic

#depth <- log10(depth)

nb <- 1
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),
     yaxt="n",xlim=c(1.69,3.5),ylim=c(0,30),yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
#lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
#lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)
axis(2,c(0,15,30),las=1,cex.axis=1.2)

nb <- 2
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)

nb <- 3
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)

###################################################################################
#### get biomass
pathname <- file.path(path, "depth_sixclasses.mat")
databio <- readMat(pathname)

forage  <- databio$BioFf
meso    <- databio$BioMf
pelpred <- databio$BioPf
batpel  <- databio$BioBaP
dem     <- databio$BioDf

depth <- as.numeric(databio$bdepth)
depth <- log10(depth)
pelprod <- as.numeric(databio$bRp)
benprod <- as.numeric(databio$bRb)
ph <- databio$param[,,1]$photic
linesize  <- 1.5

nb <- 1
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),
     yaxt="n",xlim=c(1.69,3.5),ylim=c(0,30),yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
#lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
#lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)
axis(2,c(0,15,30),las=1,cex.axis=1.2)

nb <- 2
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)

nb <- 3
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)

###################################################################################
#### get biomass
pathname <- file.path(path, "depth_nineclasses.mat")
databio <- readMat(pathname)

forage  <- databio$BioFf
meso    <- databio$BioMf
pelpred <- databio$BioPf
batpel  <- databio$BioBaP
dem     <- databio$BioDf

depth <- as.numeric(databio$bdepth)
depth <- log10(depth)
pelprod <- as.numeric(databio$bRp)
benprod <- as.numeric(databio$bRb)
ph <- databio$param[,,1]$photic
linesize  <- 1.5

nb <- 1
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),
     yaxt="n",xlim=c(1.69,3.5),ylim=c(0,30),yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
#lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
#lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)
axis(2,c(0,15,30),las=1,cex.axis=1.2)

nb <- 2
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)

nb <- 3
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("","",""),cex.axis=1)

###################################################################################
#### get biomass
pathname <- file.path(path, "depth_12classes.mat")
databio <- readMat(pathname)

forage  <- databio$BioFf
meso    <- databio$BioMf
pelpred <- databio$BioPf
batpel  <- databio$BioBaP
dem     <- databio$BioDf

depth <- as.numeric(databio$bdepth)
depth <- log10(depth)
pelprod <- as.numeric(databio$bRp)
benprod <- as.numeric(databio$bRb)
ph <- databio$param[,,1]$photic
linesize  <- 1.5

nb <- 1
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),
     yaxt="n",xlim=c(1.69,3.5),ylim=c(0,30),yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
#lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
#lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("50","250","2000"),cex.axis=1.2)
axis(2,c(0,15,30),las=1,cex.axis=1.2)

nb <- 2
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("50","250","2000"),cex.axis=1.2)

nb <- 3
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,30),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=0.75,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#8f0505",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#5c88c5")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000)),c("50","250","2000"),cex.axis=1.2)

title(xlab = "Seafloor depth (m)",
      ylab = TeX("Biomass g WW $m^{-2}$"),
      outer = TRUE, line = 2,cex.lab=1.5)
par(op)

dev.off()

#plot.new()
#legend("topleft", legend=c("large pelagic","mid-water predator","demersal","forage fish","mesopelagic fish"),
#       col=c("#8f0505", "#5c88c5", "black","#8f0505", "#5c88c5"), 
#       lwd=c(linesize,linesize,linesize,linesize,linesize), lty= c(1,1,1,2,2), cex=0.8,box.lty=0)
#save 6.90x5.27

