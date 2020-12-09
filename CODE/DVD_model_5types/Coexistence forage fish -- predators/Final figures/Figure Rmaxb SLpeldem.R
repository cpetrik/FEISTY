
rm(list=ls())
library(latex2exp)

#### script to get matlab biomass output 
library(R.matlab)
path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Coexistence forage fish -- predators/"

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Coexistence forage fish -- predators/Final figures/")
pdf("Coexistence_threespec_equalQ10.pdf",width=7.2,height=3)

op <- par(mfrow = c(1,5),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
#### get biomass
pathname <- file.path(path, "scenario_3spec.mat")
databio <- readMat(pathname)
forage  <- databio$BioFf
pelpred <- databio$BioPf
dem    <- databio$BioDem
prod <- as.numeric(databio$bRb)*10
xmax <- 50*10
linesize <- 1.5
nb  <- 1

plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
lines(dem[,nb]~prod,type="l",col="black",lwd=linesize,lty=1)
axis(1,c(0,250,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,cex.axis=1)
legend(-30,39.5,c("(a)"),  bty='n' )

nb <- 2
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
lines(dem[,nb]~prod,type="l",col="black",lwd=linesize,lty=1)
axis(1,c(0,250,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(b)"),  bty='n' )

nb <- 3
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
lines(dem[,nb]~prod,type="l",col="black",lwd=linesize,lty=1)
axis(1,c(0,250,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(c)"),  bty='n' )

nb <- 4
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
lines(dem[,nb]~prod,type="l",col="black",lwd=linesize,lty=1)
axis(1,c(0,250,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(d)"),  bty='n' )

nb <- 5
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
lines(dem[,nb]~prod,type="l",col="black",lwd=linesize,lty=1)
axis(1,c(0,250,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(e)"),  bty='n' )
legend(20,35, legend=c("epipelagic fish","large pelagics","demersal"),
       col=c("#8f0505", "#8f0505","black"), 
       lty=c(2,1,1), lwd=c(linesize,linesize,linesize),   bty='n' ,
       inset=.1,y.intersp=1, x.intersp=0.4)

title(xlab = TeX("Particle export flux (g WW $m^{-2}$ $y^{-1}$)"),
      ylab = TeX("Biomass g WW $m^{-2}$"),
      outer = TRUE, line = 2)
par(op)

dev.off()