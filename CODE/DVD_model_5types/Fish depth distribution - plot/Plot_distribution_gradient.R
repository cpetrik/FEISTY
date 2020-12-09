   # plot describing fish position in water column
    
    ## regions shallower than euphotic zone
      bottom_sh <- c(100:250)
      surf_sh <- rep(0,length(bottom_sh))
      bot_sh  <- bottom_sh
    
      # get negative
      surf_sh <- -(surf_sh)
      bot_sh <- -(bot_sh)
    
    ## regions deeper than 250 with mesopelagic / mid-water predators
      bottom_dp <- c(250:3000)
      surf_dp <- rep(0,length(bottom_dp))
      bot_dp  <- bottom_dp
      meso_dp <- rep(750,length(bottom_dp))
      meso_dp <- ifelse(meso_dp > bottom_dp, bottom_dp, meso_dp)
      dem_dp1 <- ifelse(bottom_dp-meso_dp < 1200,meso_dp,NA)
      dem_dp1 <- ifelse(bottom_dp-meso_dp >= 1200 & bottom_dp-meso_dp < 1400, bottom_dp-1200, dem_dp1)
      dem_dp2 <- ifelse(is.na(dem_dp1),bottom_dp,NA) 
      
      # get negative
      surf_dp <- -(surf_dp)
      bot_dp <- -(bot_dp)
      meso_dp <- -(meso_dp)
      dem_dp1 <- -(dem_dp1)
      dem_dp2 <- -(dem_dp2)
    
    ## set common parameters
      colb <- "blue"
      size <- 4
    
      x <- c(50:3100,50)
      x <-log10(x)
      y <- c(-50:-3100,-3100)
      
      bottom_sh <- log10(bottom_sh)
      bottom_dp <- log10(bottom_dp)
      
    # save 5 x  5.5 inch
    ## plot during the day
    
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Fish depth distribution - plot/")
      
  # plot
  pdf("Day distribution.pdf",width=5,height=5.5)   
      
    par(mfrow=c(4,3), mai = c(0.1, 0.1, 0.1, 0.1))
    
    plot(0, col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,col=colb, lwd=size)
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    lines(meso_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,col=colb, lwd=size)
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    lines(meso_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(bot_sh~bottom_sh,col=colb, lwd=size)
    lines(bot_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    # epipelagic /large pelagic
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,col=colb, lwd=size)
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,col=colb, lwd=size)
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,col=colb, lwd=size)
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    lines(meso_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    # mesopelagic/midwater
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(meso_dp~bottom_dp,col=colb, lwd=size)
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(meso_dp~bottom_dp,col=colb, lwd=size)
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(meso_dp~bottom_dp,col=colb, lwd=size)
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    # demersal
    plot(0,  col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,col=colb, lwd=size)
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0,  col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(bot_sh~bottom_sh,col=colb, lwd=size)
    lines(bot_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0,  col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(bot_sh~bottom_sh,lwd=size,col="blue")
    lines(dem_dp1~bottom_dp,lwd=size,col="blue")
    lines(dem_dp2~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    dev.off()
    
    
# now during the night

    pdf("Night distribution.pdf",width=5,height=5.5)   
    
    par(mfrow=c(4,3), mai = c(0.1, 0.1, 0.1, 0.1))
    
    # zoop s
    plot(0,  col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    # zoop L
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    #benth 
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(bot_sh~bottom_sh,lwd=size,col="blue")
    lines(bot_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    # pelagics
    plot(0,  col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    # mesopelagics/mid-water
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(meso_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0,  col="white", ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(bot_sh~bottom_sh,lwd=size,col="blue")
    lines(bot_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    plot(0, col="white",  ylim=c(-3000,0),xlim=c(2,3.5),xaxt="n", 
         yaxt="n",xlab = "",ylab = "",xaxs="i")
    polygon(x, y,border = "#d3bc5f",col="#d3bc5f")
    lines(bot_sh~bottom_sh,lwd=size,col="blue")
    lines(surf_sh~bottom_sh,lwd=size,col="blue")
    lines(bot_dp~bottom_dp,lwd=size,col="blue")
    box()
    axis(1,las=1, c(2,2.4,3,3.5),c("","","",""))
    axis(2,las=1, c(0),c(""))
    
    dev.off()
