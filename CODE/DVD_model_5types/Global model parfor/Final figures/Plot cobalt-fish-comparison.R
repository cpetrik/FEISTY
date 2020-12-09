

# model comparison of zooplankton biomass and flux to predators
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor")
  dataglob <- read.csv("datglob.csv",header=F)
  colnames(dataglob) <- c("BioZs","BioZl","BioBs","BioBl", "BioFf","BioMf","BioPf","BioBaP","BioDf")

# load input parameters
  param <- read.csv("input_parameters.csv",header=T)
  outp <- cbind(param,dataglob)
  
  setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Global model parfor/Final figures/")
  jpeg(file = "Figure_cobalt_fish_comparison.jpeg", width=7, height=5.5,units ='in', res = 500)
  op <- par(mfrow = c(1,2),
            oma = c(5,4,1,1) + 0.1,
            mar = c(3,4,1,1) + 0.1)
  
  ### plot zoopL to fish loss versus cobalt
  plot(outp$lz_bio,outp$BioZl,xlim=c(0,26),ylim=c(0,26),xlab="Large zooplankton biomass \n (COBALT)",
       ylab="Large zooplankton biomass (Fish model)",xaxt="n",yaxt="n")
  lines(c(0,25),c(0,25),col="red",lwd=2)
  axis(1,c(0,13,26))
  axis(2,c(0,13,26),las=1)
  
  zprodmod <- 1*(outp$RmaxL-outp$BioZl)
  plot(outp$RmaxL,zprodmod,ylab="Flux to fish (Fish model)",xlab="Flux to predators \n (COBALT)",
       xlim=c(0,600),ylim=c(0,600),xaxt="n",yaxt="n")
  lines(c(0,600),c(0,600),col="red",lwd=2)
  axis(1,c(0,300,600))
  axis(2,c(0,300,600),las=1)

  par(op)
  
  dev.off()