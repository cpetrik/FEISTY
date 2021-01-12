
rm(list=ls())
library(latex2exp)

#### script to get matlab biomass output 

library(R.matlab)
path <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/DVD_model_5types/Coexistence forage fish -- predators/"

setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/DVD_model_5types/Coexistence forage fish -- predators/Final figures/")
pdf("Coexistence_main_Q10unequal_60both.pdf",width=7.2,height=3)

op <- par(mfrow = c(1,5),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

#### get biomass
pathname <- file.path(path, "scenario_2spec_unequal_60both.mat")
databio <- readMat(pathname)
forage  <- databio$BioFf
pelpred <- databio$BioPf
prod <- as.numeric(databio$bRp)
xmax <- 150
linesize <- 1.5
nb  <- 1

plot(forage[,nb]~prod,type="l",lwd=linesize, col="white",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
lines(forage[,nb]~prod,type="l",col="#8f0505",lwd=linesize,lty=2)
axis(1,c(0,75,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,cex.axis=1)
legend(-15,39.5,c("(a)"),  bty='n' )

nb <- 2
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
axis(1,c(0,75,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(b)"),  bty='n' )

nb <- 3
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
axis(1,c(0,75,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(c)"),  bty='n' )

nb <- 4
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
axis(1,c(0,75,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(d)"),  bty='n' )

nb <- 5
plot(forage[,nb]~prod,type="l",lwd=linesize, col="#8f0505",
     xlab="Zooplankton production",ylab=TeX("Biomass g $m^{-2}$"),
     yaxt="n",ylim=c(0,40),xlim=c(0,xmax),lty=2,
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
lines(pelpred[,nb]~prod,type="l", col="#8f0505",lwd=linesize)
axis(1,c(0,75,xmax),cex.axis=1)
axis(2,c(0,20,40),las=1,c("","",""),cex.axis=1)
legend(-15,39.5,c("(e)"),  bty='n' )
legend(15,35, legend=c("epipelagic fish","large pelagics"),
       col=c("#8f0505", "#8f0505"), 
       lty=c(2,1), lwd=c(linesize,linesize),   bty='n' ,
       inset=.1,y.intersp=1, x.intersp=0.4)

title(xlab = TeX("Max. zooplankton prod. (g WW $m^{-2}$ $y^{-1}$)"),
      ylab = TeX("Biomass g WW $m^{-2}$"),
      outer = TRUE, line = 2)
par(op)

dev.off()

