
library(RCurl)
library(XML)
library(plyr)

setwd("C:/Users/pdvd/Online for git/Fish_foodwebs/Foodweb_plots")
dat <- read.csv("Species_groups_EwE.csv",header=T,sep=";")

# North Sea
#To Obtain a model - output
  h=basicTextGatherer()
  mymodel<- 457     
  curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  input1<-xpathSApply(data,'//group',function(x) xmlToList(x))
  nam <- unlist(input1[1,])
  id  <- unlist(input1[2,])
  nam <- data.frame(nam,as.numeric(id))
  
  feed_mat <- matrix(data=NA,nrow = nrow(nam),ncol=nrow(nam))
  col <- "diet_descr.diet.prey_seq"

  for(p in 1:nrow(nam)){
    tt <- unlist(input1[,p])
    indexcol <- which(names(tt) %in% col) 
    prey <- as.numeric(tt[indexcol])
    int <- as.numeric(tt[indexcol+1])
    feed_mat[p,prey] <- int
  }
  
  h=basicTextGatherer()
  curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client_output.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  output1<-xpathSApply(data,'//group',function(x) xmlToList(x))
  nam$Biomass <- as.numeric(unlist(output1[4,]))
  nam$Biomass[nam$Biomass < 0] <-0
  nam <- nam[order(nam$as.numeric.id.),]
  flux <- nam$Biomass %*% t(nam$Biomass)
  flux <- flux*feed_mat
  nam <- cbind(nam, dat[match(nam$nam,dat$animals), c(3,4,5)])

  agg_fin <- aggregate(nam$Biomass, list(nam$Grouping), sum)
  agg_flux <- as.data.frame(matrix(data=NA,nrow = 8, ncol=2))
  colnames(agg_flux) <- colnames(agg_fin)
  
  pred <- c("Pelagic","Pelagic","Pelagic","Pelagic","Demersal","Demersal","Demersal","Demersal")
  prey <- c("Zooplankton","Benthos","Demersal","Pelagic","Zooplankton","Benthos","Demersal","Pelagic")
  
  for (j in 1:8){
    agg_flux[j,1]  <- paste(pred[j],prey[j],sep="_")  
    subpred        <- which(nam$Grouping %in% pred[j]) 
    subprey        <- which(nam$Grouping %in% prey[j])
    agg_flux[j,2]  <- sum(flux[subpred,subprey],na.rm = T)
  }
  
  agg_fin_northsea <- rbind(agg_fin,agg_flux)

  # Peru
  #To Obtain a model - output
  h=basicTextGatherer()
  mymodel<- 658     
  curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  input1<-xpathSApply(data,'//group',function(x) xmlToList(x))
  nam <- unlist(input1[1,])
  id  <- unlist(input1[2,])
  nam <- data.frame(nam,as.numeric(id))
  
  feed_mat <- matrix(data=NA,nrow = nrow(nam),ncol=nrow(nam))
  col <- "diet_descr.diet.prey_seq"
  
  for(p in 1:nrow(nam)){
    tt <- unlist(input1[,p])
    indexcol <- which(names(tt) %in% col) 
    prey <- as.numeric(tt[indexcol])
    int <- as.numeric(tt[indexcol+1])
    feed_mat[p,prey] <- int
  }
  
  h=basicTextGatherer()
  curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client_output.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  output1<-xpathSApply(data,'//group',function(x) xmlToList(x))
  nam$Biomass <- as.numeric(unlist(output1[4,]))
  nam$Biomass[nam$Biomass < 0] <-0
  nam <- nam[order(nam$as.numeric.id.),]
  flux <- nam$Biomass %*% t(nam$Biomass)
  flux <- flux*feed_mat
  nam <- cbind(nam, dat[match(nam$nam,dat$animals), c(3,4,5)])
  
  agg_fin <- aggregate(nam$Biomass, list(nam$Grouping), sum)
  agg_flux <- as.data.frame(matrix(data=NA,nrow = 8, ncol=2))
  colnames(agg_flux) <- colnames(agg_fin)
  
  pred <- c("Pelagic","Pelagic","Pelagic","Pelagic","Demersal","Demersal","Demersal","Demersal")
  prey <- c("Zooplankton","Benthos","Demersal","Pelagic","Zooplankton","Benthos","Demersal","Pelagic")
  
  for (j in 1:8){
    agg_flux[j,1]  <- paste(pred[j],prey[j],sep="_")  
    subpred        <- which(nam$Grouping %in% pred[j]) 
    subprey        <- which(nam$Grouping %in% prey[j])
    agg_flux[j,2]  <- sum(flux[subpred,subprey],na.rm = T)
  }
  
  agg_fin_Peru <- rbind(agg_fin,agg_flux)
  
  # Scot
  #To Obtain a model - output
  h=basicTextGatherer()
  mymodel<- 443     
  curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  input1<-xpathSApply(data,'//group',function(x) xmlToList(x))
  nam <- unlist(input1[1,])
  id  <- unlist(input1[2,])
  nam <- data.frame(nam,as.numeric(id))
  
  feed_mat <- matrix(data=NA,nrow = nrow(nam),ncol=nrow(nam))
  col <- "diet_descr.diet.prey_seq"
  
  for(p in 1:nrow(nam)){
    tt <- unlist(input1[,p])
    indexcol <- which(names(tt) %in% col) 
    prey <- as.numeric(tt[indexcol])
    int <- as.numeric(tt[indexcol+1])
    feed_mat[p,prey] <- int
  }
  
  h=basicTextGatherer()
  curlPerform(url = paste('http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client_output.php?no_model=',mymodel,sep=''),writefunction=h$update,verbose=TRUE)
  data<-xmlTreeParse(h$value(),useInternalNodes=TRUE)
  output1<-xpathSApply(data,'//group',function(x) xmlToList(x))
  nam$Biomass <- as.numeric(unlist(output1[4,]))
  nam$Biomass[nam$Biomass < 0] <-0
  nam <- nam[order(nam$as.numeric.id.),]
  flux <- nam$Biomass %*% t(nam$Biomass)
  flux <- flux*feed_mat
  nam <- cbind(nam, dat[match(nam$nam,dat$animals), c(3,4,5)])
  
  agg_fin <- aggregate(nam$Biomass, list(nam$Grouping), sum)
  agg_flux <- as.data.frame(matrix(data=NA,nrow =15, ncol=2))
  colnames(agg_flux) <- colnames(agg_fin)
  
  pred <- c("Pelagic","Pelagic","Pelagic","Pelagic","Pelagic","Demersal","Demersal","Demersal","Demersal","Demersal",
            "Midwater","Midwater","Midwater","Midwater","Midwater")
  prey <- c("Zooplankton","Benthos","Demersal","Pelagic","Midwater","Zooplankton","Benthos","Demersal","Pelagic","Midwater",
            "Zooplankton","Benthos","Demersal","Pelagic","Midwater")
  
  for (j in 1:15){
    agg_flux[j,1]  <- paste(pred[j],prey[j],sep="_")  
    subpred        <- which(nam$Grouping %in% pred[j]) 
    subprey        <- which(nam$Grouping %in% prey[j])
    agg_flux[j,2]  <- sum(flux[subpred,subprey],na.rm = T)
  }
  
  agg_fin_Scot <- rbind(agg_fin,agg_flux)
  
  rm(list=setdiff(ls(), c("agg_fin_northsea","agg_fin_Peru","agg_fin_Scot","path")))
  
  
  