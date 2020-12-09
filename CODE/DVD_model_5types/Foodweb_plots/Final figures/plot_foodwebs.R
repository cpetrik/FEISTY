  
# save 7.2 x 6.5 portrait

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/Final figures/")

# plot on one page
pdf("Fish foodweb plots.pdf",width=6.5,height=6) 

op <- par(mfrow = c(3,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(1,3.5,1.5,1) + 0.1)

# inkscape unable to load circle size properly with par --> redo 

#pdf("Fish foodweb plots_sep.pdf",width=5,height=4) 

############figure A food web shallow temperate ######################
  rm(list=ls())
  
  library(latex2exp)
  library(R.matlab)

  path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/" 

 # get info
  pathname <- file.path(path, "region_1.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  position <- as.numeric( -(foodweb$param[,,1]$avlocDay + foodweb$param[,,1]$avlocNight)/2)

  # check amount of biomass and onlygroups with biomass total > 0.1 are included
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])

  # define parameters 
  idx <- which(biomass > 0.1)
  mat <- matrix(data=0,nrow=30)
  mat[idx,1]<- 1
  biomass <- log10(biomass+1)*20
  bio_cex <- sqrt(((4*pi*1^2)*biomass)/ (4*pi))
  nointer <- 10^-2
  breakinter <- c(0,0.1,1,100)
  flux <- as.matrix(foodweb$mortpr)
  linesize <- c(1,1,2)
  linetype <- c(3, 1, 1)
  pos <- c(position[1:2],position[3:4]-2,position[5:8]+0,position[9:24]-1,
           position[25:26]-2 , position[27:30]+0)
  
  # create plot
  plot(log10(upsize[idx]),position[idx],cex=bio_cex[idx],ylim=c(-55,5),yaxt="n",
       xlim=c(-4,5),col="white", xlab = "",xaxt="n",ylab="",main="shallow shelf")
  axis(1,c(-4,-2,0,2,4),labels = c("","","","",""))
  axis(2,c(0,-50),labels=c(0,50),las=2,cex.axis=1.4)
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[13:18]),pos[13:18],col="#fa2306",cex=bio_cex[13:18],pch=16)
  points(log10(upsize[25:30]),pos[25:30],col="black",cex=bio_cex[25:30],pch=16)

  # add interactions for the different groups
  for (nb in 25:30){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac >nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="black")
  }}}

  
  for (nb in 13:18){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#fa2306")
      }}}
  
  for (nb in 5:8){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#8f0505")
      }}}

  nb <- 1
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer){
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                              labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#2ca25f")
  }}
    
  nb <- 2
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#99d8c9")
  }}

  nb <- 3
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#d95f0e")
  }}

  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[13:18]),pos[13:18],col="#fa2306",cex=bio_cex[13:18],pch=16)
  points(log10(upsize[25:30]),pos[25:30],col="black",cex=bio_cex[25:30],pch=16)
  box()
  
  ##### figure 2
  rm(list=ls())
  
  library(latex2exp)
  library(R.matlab)
  path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/"
  
  # get info
  pathname <- file.path(path, "region_2.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  position <- as.numeric( -(foodweb$param[,,1]$avlocDay + foodweb$param[,,1]$avlocNight)/2)
  
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])
  
  # define parameters 
  idx <- which(biomass > 0)
  mat <- matrix(data=0,nrow=30)
  mat[idx,1]<- 1
  biomass <- log10(biomass+1)*20
  bio_cex <- sqrt(((4*pi*1^2)*biomass)/ (4*pi))
  nointer <- 10^-2
  breakinter <- c(0,0.1,1,100)
  flux <- as.matrix(foodweb$mortpr)
  linesize <- c(1,1,2)
  linetype <- c(3, 1, 1)
  pos <- c(position[1],position[2],position[3:4]+0,position[5:8]+0,position[9:12]-25,
           position[13:16]-10,position[17:18],position[19:20],position[21:22],
           position[23:24], position[25:26]-25, position[27:30])

  # start plotting
  plot(log10(upsize[idx]),position[idx],cex=bio_cex[idx],ylim=c(-800,100),yaxt="n",
       xlim=c(-4,5),col="white", xlab = "",xaxt="n",ylab="",main="slope")
  axis(1,c(-4,-2,0,2,4),labels = c("","","","",""))
  axis(2,c(0,-750),labels=c(0,750),las=2,cex.axis=1.4)           
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  #points(log10(upsize[4]),pos[4],col="#fec44f",cex=bio_cex[4],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)
  points(log10(upsize[13:18]),pos[13:18],col="#fa2306",cex=bio_cex[13:18],pch=16)
  points(log10(upsize[19:24]),pos[19:24],col="#bcbddc",cex=bio_cex[19:24],pch=16)
  points(log10(upsize[25:30]),pos[25:30],col="#000000",cex=bio_cex[25:30],pch=16)

  # add interactions
  nb <- 1
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#2ca25f")
  }}

  nb <- 2
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer ) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#99d8c9")
  }}

  nb <- 3
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#d95f0e")
  }}

  for (nb in 5:8){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = c(0,0.01,1,100), 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#8f0505")
      }}}
  

  for (nb in 9:12){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#5c88c5")
    }}}

  for (nb in 13:18){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#fa2306")
      }}}
  
  
  for (nb in 19:24){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#bcbddc")
    }}}
  
  for (nb in 25:30){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#000000")
    }}}

  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)
  points(log10(upsize[13:18]),pos[13:18],col="#fa2306",cex=bio_cex[13:18],pch=16)
  points(log10(upsize[19:24]),pos[19:24],col="#bcbddc",cex=bio_cex[19:24],pch=16)
  points(log10(upsize[25:30]),pos[25:30],col="#000000",cex=bio_cex[25:30],pch=16)
  box()
  
#### figure 3
  rm(list=ls())
  
  library(latex2exp)
  library(R.matlab)
  library(plotrix)
  path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/"
  
  #### food web deep temperate
  # get info
  pathname <- file.path(path, "region_3.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  position <- as.numeric( -(foodweb$param[,,1]$avlocDay + foodweb$param[,,1]$avlocNight)/2)
  
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])
  
  # define parameters 
  idx <- which(biomass > 0)
  mat <- matrix(data=0,nrow=30)
  mat[idx,1]<- 1
  biomass <- log10(biomass+1)*20
  bio_cex <- sqrt(((4*pi*1^2)*biomass)/ (4*pi))
  nointer <- 10^-2
  breakinter <- c(0,0.1,1,100)
  flux <- as.matrix(foodweb$mortpr)
  linesize <- c(1,1,2)
  linetype <- c(3, 1, 1)
  pos <- c(position[1],position[2],position[3:4]+2400,position[5:8]+0,position[9:12]-25,
           position[13:16]-10,position[17:18],position[19:20],position[21:22],
           position[23:24], position[25:26]-25, position[27:30])
  
  # start plotting
  plot(log10(upsize[idx]),position[idx],cex=bio_cex[idx],ylim=c(-1500,100),yaxt="n",
       xlim=c(-4,5),col="white", xlab = "",xaxt="n",ylab="",main="upwelling")
  axis(1,c(-4,-2,0,2,4),labels = c("","","","",""))
  axis(2,c(0,-700,-1400),labels=c(0,700,3800),las=2,cex.axis=1.4)           
  axis.break(2, -1200, style = "zigzag")
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)
  points(log10(upsize[13:18]),pos[13:18],col="#fa2306",cex=bio_cex[13:18],pch=16)
  points(log10(upsize[19:24]),pos[19:24],col="#bcbddc",cex=bio_cex[19:24],pch=16)
  points(log10(upsize[25:30]),pos[25:30],col="#000000",cex=bio_cex[25:30],pch=16)

  # add interactions
  nb <- 1
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#2ca25f")
    }}
  
  nb <- 2
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer ) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#99d8c9")
    }}
  
  nb <- 3
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#d95f0e")
    }}
  

  for (nb in 5:8){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#8f0505")
      }}}
  
  
  for (nb in 9:12){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#5c88c5")
      }}}
  
  for (nb in 13:18){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#fa2306")
      }}}
  
  
  for (nb in 19:24){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#bcbddc")
      }}}
  
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)
  points(log10(upsize[13:18]),pos[13:18],col="#fa2306",cex=bio_cex[13:18],pch=16)
  points(log10(upsize[19:24]),pos[19:24],col="#bcbddc",cex=bio_cex[19:24],pch=16)
  points(log10(upsize[25:30]),pos[25:30],col="#000000",cex=bio_cex[25:30],pch=16)
  box()
  axis.break(2, -1200, style = "zigzag")

#### figure 4 
  rm(list=ls())
  
  library(latex2exp)
  library(R.matlab)
  library(plotrix)
  path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/"
  
  #### oligo 1
  # get info
  pathname <- file.path(path, "region_4.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  position <- as.numeric( -(foodweb$param[,,1]$avlocDay + foodweb$param[,,1]$avlocNight)/2)
  
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])
  
  # define parameters 
  idx <- which(biomass > 0)
  mat <- matrix(data=0,nrow=30)
  mat[idx,1]<- 1
  biomass <- log10(biomass+1)*20
  bio_cex <- sqrt(((4*pi*1^2)*biomass)/ (4*pi))
  nointer <- 10^-2
  breakinter <- c(0,0.1,1,100)
  flux <- as.matrix(foodweb$mortpr)
  linesize <- c(1,1,2)
  linetype <- c(3, 1, 1)
  pos <- c(position[1],position[2],position[3:4]+2400,position[5:8]+0,position[9:12]-25,
           position[13:16]-10,position[17:18],position[19:20],position[21:22],
           position[23:24], position[25:26]-25, position[27:30])
  
  # start plotting
  plot(log10(upsize[idx]),position[idx],cex=bio_cex[idx],ylim=c(-1500,100),yaxt="n",
       xlim=c(-4,5),col="white", xlab = "",xaxt="n",ylab="",main="oligo 1")
  axis(1,c(-4,-2,0,2,4),labels = c("","","","",""))
  axis(2,c(0,-700,-1400),labels=c(0,700,3800),las=2,cex.axis=1.4)           
  axis.break(2, -1200, style = "zigzag")
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)

  # add interactions
  nb <- 1
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#2ca25f")
    }}
  
  nb <- 2
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer ) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#99d8c9")
    }}
  
  nb <- 3
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#d95f0e")
    }}
  
  for (nb in 5:8){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#8f0505")
      }}}
  
  
  for (nb in 9:12){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#5c88c5")
      }}}

  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[5:8]),pos[5:8],col="#8f0505",cex=bio_cex[5:8],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)

  box()
  axis.break(2, -1200, style = "zigzag")
  
  #### figure 5 
  rm(list=ls())
  
  library(latex2exp)
  library(R.matlab)
  library(plotrix)
  path <- "C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots/"
  
  #### oligo 2
  # get info
  pathname <- file.path(path, "region_5.mat")
  foodweb <- readMat(pathname)
  upsize <- as.numeric(foodweb$param[,,1]$wc)
  biomass  <- as.numeric(foodweb$yend)
  position <- as.numeric( -(foodweb$param[,,1]$avlocDay + foodweb$param[,,1]$avlocNight)/2)
  
  spel <- sum(biomass[5:8])
  mpel <- sum(biomass[9:12])
  lpel <- sum(biomass[13:18])
  bpel <- sum(biomass[19:24])
  dem <- sum(biomass[25:30])
  
  # define parameters 
  idx <- which(biomass > 0)
  mat <- matrix(data=0,nrow=30)
  mat[idx,1]<- 1
  biomass <- log10(biomass+1)*20
  bio_cex <- sqrt(((4*pi*1^2)*biomass)/ (4*pi))
  nointer <- 10^-2
  breakinter <- c(0,0.1,1,100)
  flux <- as.matrix(foodweb$mortpr)
  linesize <- c(1,1,2)
  linetype <- c(3, 1, 1)
  pos <- c(position[1],position[2],position[3:4]+3800,position[5:8]+0,position[9:12]-25,
           position[13:16]-10,position[17:18],position[19:20],position[21:22],
           position[23:24], position[25:26]-25, position[27:30])
  
  # start plotting
  plot(log10(upsize[idx]),position[idx],cex=bio_cex[idx],ylim=c(-1500,100),yaxt="n",
       xlim=c(-4,5),col="white", xlab = "",xaxt="n",ylab="",main="oligo 2")
  axis(1,c(-4,-2,0,2,4),labels = c(TeX("$10^{-4}$"),0.01,1,100,10000),cex.axis=1.4)
  axis(2,c(0,-700,-1400),labels=c(0,700,5200),las=2,cex.axis=1.4)           
  axis.break(2, -1200, style = "zigzag")
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)

  # add interactions
  nb <- 1
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#2ca25f")
    }}
  
  nb <- 2
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer ) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#99d8c9")
    }}
  
  nb <- 3
  for(i in 1:30){
    interac <- flux[i,nb]*mat[i]
    if (interac > nointer) {
      lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                        labels =c(1:3))))
      lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
            lty=linetype[lz],col="#d95f0e")
    }}
  
 
  
  for (nb in 9:12){
    for(i in 1:30){
      interac <- flux[i,nb]*mat[i]
      if (interac > nointer) {
        lz <- as.numeric(as.character(cut(interac, breaks = breakinter, 
                                          labels =c(1:3))))
        lines(log10(upsize[c(nb,i)]),pos[c(nb,i)],lwd=linesize[lz],
              lty=linetype[lz],col="#5c88c5")
      }}}
  
  points(log10(upsize[1]),pos[1],col="#2ca25f",cex=bio_cex[1],pch=16)
  points(log10(upsize[2]),pos[2],col="#99d8c9",cex=bio_cex[2],pch=16)
  points(log10(upsize[3]),pos[3],col="#d95f0e",cex=bio_cex[3],pch=16)
  points(log10(upsize[9:12]),pos[9:12],col="#5c88c5",cex=bio_cex[9:12],pch=16)
  
  box()
  axis.break(2, -1200, style = "zigzag")
  
  title(xlab = "Central weight (g WW)",
        ylab = "Depth (m)",
        outer = TRUE, line = 2,cex.lab=1.4)
  
  dev.off()
