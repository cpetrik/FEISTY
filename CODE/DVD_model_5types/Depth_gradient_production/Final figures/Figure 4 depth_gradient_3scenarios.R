
rm(list=ls())
library(latex2exp)

#### script to get matlab biomass output 

library(R.matlab)
path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Depth_gradient_production"

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
linesize  <- 3

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Depth_gradient_production/Final figures")
pdf("Figure_4.pdf",width=3.3,height=5.5)

op <- par(mfrow = c(3,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

nb <- 3
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,22),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#fa2306",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#bcbddc")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000),log10(3000)),c("","","",""),cex.axis=1)
axis(2,c(0,11,22),las=1,cex.axis=1.2)
legend(log10(1700),22,c("(a)"),  bty='n' )

nb <- 2
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth",ylab=TeX("Biomass gr $m^{-2}$"),yaxt="n",ylim=c(0,14),
     yaxs="i",xaxt="n",xaxs="i",cex.lab=1,xlim=c(1.69,3.5))
lines(pelpred[,nb]~depth,type="l", col="#fa2306",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#bcbddc")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000),log10(3000)),c("","","",""),cex.axis=1)
axis(2,c(0,7,14),las=1,cex.axis=1.2)
legend(log10(1700),14,c("(b)"),  bty='n' )

nb <- 1
plot(forage[,nb]~depth,type="l",lwd=linesize, col="#8f0505",lty=2,
     xlab="seafloor depth (m)",ylab=TeX("Biomass gr $m^{-2}$"),
     yaxt="n",xlim=c(1.69,3.5),ylim=c(0,6),yaxs="i",xaxt="n",xaxs="i",cex.lab=1)
#lines(pelpred[,nb]~depth,type="l", col="#fa2306",lwd=linesize)
lines(meso[,nb]~depth,type="l",lwd=linesize, col="#5c88c5",lty=2)
#lines(batpel[,nb]~depth,type="l",lwd=linesize, col="#bcbddc")
lines(dem[,nb]~depth,type="l",lwd=linesize, col="black")
axis(1,c(log10(50),log10(250),log10(2000),log10(3000)),c("50","250","2000",""),cex.axis=1.2)
axis(2,c(0,3,6),las=1,cex.axis=1.2)
legend(log10(1700),5.9,c("(c)"),  bty='n' )


title(xlab = "Seafloor depth (m)",
      ylab = TeX("Biomass g WW $m^{-2}$"),
      outer = TRUE, line = 2,cex.lab = 1.5)
par(op)

# save 3.15 x 5

dev.off()

# detrital change

pdf("Figure_4_detrital change.pdf",width=3.7,height=2.4)
op <- par(mfrow = c(1,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

tt <- c(50:3000)

y2 <- (tt/150)^-0.86
y2[y2>1] <- 1

plot(y2~log10(tt), ylim =c(0,1),xlim=c(1.69,3.5),yaxt="n",xaxt="n",xaxs="i",yaxs="i",type="l",lwd=2,col="white",bty="n")
polygon(c(1,y2,0,0,1)~c(1.69,log10(tt),3.477,1.69,1.69),col="#d3bc5f")

lines(c(log10(150),log10(150)),c(0,1),lty=3)
axis(1,c(log10(50),log10(250),log10(2000),log10(3000)),c("","","",""),cex.axis=1)
axis(2,c(0,1),las=1)

dev.off()



