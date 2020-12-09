rm(list=ls())

#### script to get medium and large zooplankton energy flux to fish and detritus flux to seabed
#### checked to use as input in the fishmodel 

### input = "mg WW m-2 day-1"
### output = "g WW m-2 y-1"

   library(R.matlab)
   library(latex2exp)

# path to data charlie
   path <- "H:/Werk/Global drivers BP ratio/Final_analysis_manuscript/Environmental_data_Ecoregions/Charlie data" # outside github

#### get higher predator loss medium zooplankton 
   pathname <- file.path(path, "hploss_100_mdz_1deg_ESM26_5yr_clim_191_195.mat")
   datamdz <- readMat(pathname)
   
   mdz_Month<-list(datamdz$hploss.100.mdz[1,,], datamdz$hploss.100.mdz[2,,],datamdz$hploss.100.mdz[3,,],datamdz$hploss.100.mdz[4,,],
                   datamdz$hploss.100.mdz[5,,], datamdz$hploss.100.mdz[6,,],datamdz$hploss.100.mdz[7,,],datamdz$hploss.100.mdz[8,,],
                   datamdz$hploss.100.mdz[9,,], datamdz$hploss.100.mdz[10,,],datamdz$hploss.100.mdz[11,,],datamdz$hploss.100.mdz[12,,])
   
   mdz_output<-Reduce('+', mdz_Month)/12
   mdz_output<-mdz_output*365/1000*9

#### get higher predator loss large zooplankton
   pathname <- file.path(path, "hploss_100_lgz_1deg_ESM26_5yr_clim_191_195.mat")
   datalgz <- readMat(pathname)
   
   lgz_Month<-list(datalgz$hploss.100.lgz[1,,], datalgz$hploss.100.lgz[2,,],datalgz$hploss.100.lgz[3,,],datalgz$hploss.100.lgz[4,,],
                   datalgz$hploss.100.lgz[5,,], datalgz$hploss.100.lgz[6,,],datalgz$hploss.100.lgz[7,,],datalgz$hploss.100.lgz[8,,],
                   datalgz$hploss.100.lgz[9,,], datalgz$hploss.100.lgz[10,,],datalgz$hploss.100.lgz[11,,],datalgz$hploss.100.lgz[12,,])
   
   lgz_output<-Reduce('+', lgz_Month)/12
   lgz_output<-lgz_output*365/1000*9

### now get biomass medium zooplankton
### get zooplankton biomass ###
   pathname <- file.path(path, "mdz_100_1deg_ESM26_5yr_clim_191_195.mat")
   datamdz <- readMat(pathname)
   
   mdzMonth<-list(datamdz$mdz.100[1,,], datamdz$mdz.100[2,,],datamdz$mdz.100[3,,], datamdz$mdz.100[4,,],
                  datamdz$mdz.100[5,,], datamdz$mdz.100[6,,],datamdz$mdz.100[7,,], datamdz$mdz.100[8,,],
                  datamdz$mdz.100[9,,], datamdz$mdz.100[10,,],datamdz$mdz.100[11,,], datamdz$mdz.100[12,,])
   
   mdz_biomass<-Reduce('+', mdzMonth)/12/1000*9
   
   
   pathname <- file.path(path, "lgz_100_1deg_ESM26_5yr_clim_191_195.mat")
   datalgz <- readMat(pathname)
   
   lgzMonth<-list(datalgz$lgz.100[1,,], datalgz$lgz.100[2,,],datalgz$lgz.100[3,,], datalgz$lgz.100[4,,],
                  datalgz$lgz.100[5,,], datalgz$lgz.100[6,,],datalgz$lgz.100[7,,], datalgz$lgz.100[8,,],
                  datalgz$lgz.100[9,,], datalgz$lgz.100[10,,],datalgz$lgz.100[11,,], datalgz$lgz.100[12,,])
   
   lgz_biomass<-Reduce('+', lgzMonth)/12/1000*9

# make plot and compare 
   # path to matlab output coexistence epipelagic and large pelagic predator
   pathmat <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Max versus realized production/"
   
   setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Max versus realized production/Final figures/")
   
   jpeg(file = "Figure_cobalt_fish_predation.jpeg", width=7, height=5.5,units ='in', res = 500)

   op <- par(mfrow = c(1,2),
             oma = c(5,4,0,0) + 0.1,
             mar = c(1,1,1,1) + 0.1)
   
   plot(mdz_output~mdz_biomass,xlab=TeX("Biomass g $m^{-2}$"), ylab=TeX("Loss to predators (gr WW $m^{-2}$ $y^{-1}$)"),
        main="small zooplankton",las=1,ylim=c(0,600),xlim=c(0,18), col="grey")
   
   legend(0,620, legend=c(TeX("F mod. $0^{o}$C"), TeX("F mod. $10^{o}$C"), TeX("F mod. $20^{o}$C"), TeX("F mod. $30^{o}$C"),"COBALT output"),
          col=c("blue","black","orange","red","grey"), 
          lty=c(2,1,3,4,NA),pch=c(NA,NA,NA,NA,1), lwd=c(2,2,2,2,2), box.lty=0 , 
          inset=.1,y.intersp=1.5, x.intersp=0.4,  bty='n')
   legend(14,620,c("(a)"),  bty='n' )

#### load matlab output from epipelagic - large pelagic coexistence
   dat_trop <- read.csv(file = paste(pathmat,"pel_prod_temp.csv",sep="/"),header=F)
   colnames(dat_trop) <- c("FF0","FF10","FF20","FF30","PF0","PF10","PF20","PF30",
                           "ZL0","ZL10","ZL20","ZL30","ZS0","ZS10","ZS20","ZS30","bRp")
   dat_trop$prodZL0 <- 1*(dat_trop$bRp - dat_trop$ZL0)
   dat_trop$prodZL10 <- 1*(dat_trop$bRp - dat_trop$ZL10)
   dat_trop$prodZL20 <- 1*(dat_trop$bRp - dat_trop$ZL20)
   dat_trop$prodZL30 <- 1*(dat_trop$bRp - dat_trop$ZL30)
   dat_trop$prodZS0 <- 1*(dat_trop$bRp - dat_trop$ZS0)
   dat_trop$prodZS10 <- 1*(dat_trop$bRp - dat_trop$ZS10)
   dat_trop$prodZS20 <- 1*(dat_trop$bRp - dat_trop$ZS20)
   dat_trop$prodZS30 <- 1*(dat_trop$bRp - dat_trop$ZS30)

   lines(dat_trop$prodZS0~dat_trop$ZS0, col="blue",pch=16,lwd=3,lty=2)
   lines(dat_trop$prodZS10~dat_trop$ZS10, col="black",pch=16,lwd=3,lty=1)
   lines(dat_trop$prodZS20~dat_trop$ZS20, col="orange",pch=16,lwd=3,lty=3)
   lines(dat_trop$prodZS30~dat_trop$ZS30, col="red",pch=16,lwd=3,lty=4)

   plot(lgz_output~lgz_biomass,xlab=TeX("Biomass g $m^{-2}$"), ylab="",
        main="large zooplankton",las=1,ylim=c(0,600),xlim=c(0,18),yaxt="n", col="grey")
   lines(dat_trop$prodZL0~dat_trop$ZL0, col="blue",pch=16,lwd=3,lty=2)
   lines(dat_trop$prodZL10~dat_trop$ZL10, col="black",pch=16,lwd=3,lty=1)
   lines(dat_trop$prodZL20~dat_trop$ZL20, col="orange",pch=16,lwd=3,lty=3)
   lines(dat_trop$prodZL30~dat_trop$ZL30, col="red",pch=16,lwd=3,lty=4)
   points(dat_trop$prodZL10[7]~dat_trop$ZL10[7], col="black",pch=16,cex=1.5)
   points(dat_trop$prodZL10[25]~dat_trop$ZL10[25], col="black",pch=16,cex=1.5)
   points(dat_trop$prodZL10[46]~dat_trop$ZL10[46], col="black",pch=16,cex=1.5)

   legend(14,620,c("(b)"),  bty='n' )

   title(ylab = TeX("Loss to predators (g WW $m^{-2}$ $y^{-1}$)"),
         xlab = TeX("Biomass g WW $m^{-2}$"),
         outer = TRUE, line = 2)
   par(op)
   
   dev.off()

# now get detrital seabed data and depth
   ##### get detritus flux at seabed
   pathname <- file.path(path, "fcdet_btm_1deg_ESM26_5yr_clim_191_195.mat")
   datafcdet <- readMat(pathname)
   
   fcdet_Month<-list(datafcdet$fcdet.btm[1,,], datafcdet$fcdet.btm[2,,],datafcdet$fcdet.btm[3,,],datafcdet$fcdet.btm[4,,],
                     datafcdet$fcdet.btm[5,,], datafcdet$fcdet.btm[6,,],datafcdet$fcdet.btm[7,,],datafcdet$fcdet.btm[8,,],
                     datafcdet$fcdet.btm[9,,], datafcdet$fcdet.btm[10,,],datafcdet$fcdet.btm[11,,],datafcdet$fcdet.btm[12,,])
   
   fcdet_output<-Reduce('+', fcdet_Month)/12
   fcdet_output<-fcdet_output*365/1000
   
   fcdet_output <- fcdet_output*9 ## carbon to ww

### load depth data Charlie
   path <- "H:/Werk/Global drivers BP ratio/Final_analysis_manuscript/Environmental_data_Ecoregions/Charlie data" # outside github
   pathname <- file.path(path, "depth_1deg_ESM26_lat_lon.mat")
   depth <- readMat(pathname)
   depth_mat <- depth$depth
   
   fcdet_shal <- fcdet_output
   fcdet_shal[depth_mat > 100] <- NA
   quantile(fcdet_shal,c(0.15,0.854),na.rm=T)
   
   lgz_shal <- lgz_output
   lgz_shal[depth_mat > 100] <- NA
   
   plot(fcdet_output,lgz_shal)
   
   Flux <- as.vector(t(fcdet_output))
   Lprod <- as.vector(t(lgz_shal))
   
   # plot detritus flux
   jpeg(file = "Figure sub_detritus_flux.jpeg", width=5, height=4.5,units ='in', res = 500)
   
   op <- par(mfrow = c(1,1),
             oma = c(5,4,0,0) + 0.1,
             mar = c(1,1,1,1) + 0.1)
   
   plot(Flux~Lprod, col="grey",las=1)
   low <- dat_trop$prodZL10[7]
   med <- dat_trop$prodZL10[25]
   high <- dat_trop$prodZL10[46]

   # some stats
   #summary(lm(Flux~Lprod))
   x <- c(0:600)
   y <- 121.35663 + 2.58560*x
   lines(y~x, col="blue",lwd=2)
   
   x <- c(low,med,high)
   y <- 121.35663 + 2.58560*x
   points(y~x, pch=16, col="black", cex=2)
   
   legend(300,400, legend=c("linear regression","COBALT output"),
          col=c("blue","grey"), 
          lty=c(1,NA),pch=c(NA,1), lwd=c(2,2), box.lty=0 , 
          inset=.1,y.intersp=1.5, x.intersp=0.4)

   title(ylab = TeX("Detrital flux to the seabed (g WW $m^{-2}$ $y^{-1}$)"),
         xlab = TeX("Loss to predators (g WW $m^{-2}$ $y^{-1}$)"),
         outer = TRUE, line = 2)
   par(op)
   dev.off()
