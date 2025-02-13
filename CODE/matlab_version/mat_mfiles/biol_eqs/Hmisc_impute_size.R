################################################################################

# Impute missing values of weight and length 
# Hmisc imputation package
# RAM v3.0 and ICES stocks species 
# ALL species, not just the ones with the time series we deemed acceptable
# Plus additional species from RAM v1.0 not represented here (lots of sharks)

################################################################################

rm(list=ls())

library(Hmisc)
library(VIM)

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/RAM_Legacy/ICES_RAMv3.0ado/")


### load data
#All species from all stocks
traits <- read.csv("ICES_RAMv3.0_ALL_spp_mean_traits.csv",sep=",", header = T, stringsAsFactors = F)
misv1 <- read.csv("RAMv1.0_extra_spp_mean_FBtraits_filled_tax_fec_mat.csv",sep=",", header = T, stringsAsFactors = F)

### Combine 
# Remove Gadus macrocephalus and Sparus aurata from "extra"
mtsc <- traits$StockCode
mvsc <- misv1$StockCode
rep <- intersect(mtsc,mvsc)
mvsc2 <- setdiff(mvsc,rep)
misv12 <- as.data.frame(misv1[(misv1$StockCode %in% mvsc2),])
mt<-names(traits)
mv<-names(misv12)
setdiff(mt,mv)
misv12$RAM_SciName <- NA
misv12$PCI <- NA
misv12 <- as.data.frame(misv12[,mt])
mtraits <- rbind(traits,misv12) #477 stocks


### Reduce to traits for imputation
# Make StockCode a factor
mtraits$StockCode <- as.factor(mtraits$StockCode)

length <- c("TLinf","TLmax","TLBoth")
weight <- c("Weight","Wmax","Winf")

col <- c("StockCode","sciname","Order",length, weight)
smtraits <- mtraits[,col]

### Missingness
aggr_plot <- aggr(smtraits, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(smtraits), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))
##Least missing:
#Weight 0.21593291
#TLBoth   0.03983229

sco <- which(smtraits$Order=="Scorpaeniformes") #50
raj <- which(smtraits$Order=="Rajiformes")      #3
squ <- which(smtraits$Order=="Squaliformes")    #9
chn <- c(raj,squ)
fish <- which(smtraits$Order!="Scorpaeniformes")
bony <- setdiff(fish,chn)                       #415


# Log transform
smtraits[,10:15] <- log(smtraits[,4:9])
names(smtraits)[10:15] <- c("logTLinf","logTLmax","logTLBoth",
                            "logWeight","logWmax","logWinf")


#------------------------ IMPUTE TRAITS by Order ----------------------------------

### Non-scorp, Non-chondr
dbony <- as.data.frame(smtraits[bony,])
bony_fec <- aregImpute(~ logTLinf+logTLmax+logTLBoth+logWeight+logWmax+logWinf, 
                       data = dbony, n.impute = 10, nk=0)
bony_fec
#R-squares for Predicting Non-Missing Values for Each Variable
#Using Last Imputations of Predictors
#logTLinf  logTLmax logTLBoth logWeight   logWmax   logWinf 
#0.959     0.943     0.969     0.925      0.790     0.948

## Take mean
# get a list of those variables with imputed values
bn.imp.vars <- names(bony_fec$imputed)
# compute mean imputed value for each variable
# and extract the original indices of the missing data, by variable
bn.imp.vars.mean <- lapply(bony_fec$imputed, function(i) apply(i, 1, mean))
bn.imp.vars.idx <- lapply(bn.imp.vars.mean, function(i) as.integer(names(i)))
# copy origial data
bon.no.na <- dbony
# loop over imputed variables
for(i in bn.imp.vars)
{
  print(i)
  # get the mean imputations for this variable
  bn.imp.i <- bn.imp.vars.mean[[i]]
  # get the original indices for NA
  bn.idx.i <- which(is.na(bon.no.na[[i]]))
  # replace original NA with imputed values
  bon.no.na[bn.idx.i, i] <- bn.imp.i
}



### Scorp ------
dsco <- as.data.frame(smtraits[sco,])
sco_fec <- aregImpute(~ logTLinf+logTLmax+logTLBoth+logWeight+logWmax+logWinf, 
                      data = dsco, n.impute = 10, nk=0)
sco_fec
#R-squares for Predicting Non-Missing Values for Each Variable
#Using Last Imputations of Predictors
#logTLinf  logTLmax logTLBoth logWeight   logWmax   logWinf 
#0.809     0.955     0.748     0.942      1.000     1.000 

### Combine imputations (take mean) 
##This takes the mean
# get a list of those variables with imputed values
sc.imp.vars <- names(sco_fec$imputed)
# compute mean imputed value for each variable
# and extract the original indices of the missing data, by variable
sc.imp.vars.mean <- lapply(sco_fec$imputed, function(i) apply(i, 1, mean))
sc.imp.vars.idx <- lapply(sc.imp.vars.mean, function(i) as.integer(names(i)))
# copy origial data
sco.no.na <- dsco
# loop over imputed variables
for(i in sc.imp.vars)
{
  print(i)
  # get the mean imputations for this variable
  sc.imp.i <- sc.imp.vars.mean[[i]]
  # get the original indices for NA
  sc.idx.i <- which(is.na(sco.no.na[[i]]))
  # replace original NA with imputed values
  sco.no.na[sc.idx.i, i] <- sc.imp.i
}




### Chondrich ------
dchon <- as.data.frame(smtraits[chn,])
chon_fec <- aregImpute(~ logTLinf+logTLmax+logTLBoth+logWeight+logWinf, 
                       data = dchon, n.impute = 10, nk=0)
chon_fec
#R-squares for Predicting Non-Missing Values for Each Variable
#Using Last Imputations of Predictors
#logTLinf logWeight   logWinf 
#1        -Inf         1


### Combine imputations (take mean) 
# get a list of those variables with imputed values
ch.imp.vars <- names(chon_fec$imputed)[c(1,4,5)]
# compute mean imputed value for each variable
# and extract the original indices of the missing data, by variable
ch.imp.vars.mean <- lapply(chon_fec$imputed[c(1,4,5)], function(i) apply(i, 1, mean))
ch.imp.vars.idx <- lapply(ch.imp.vars.mean, function(i) as.integer(names(i)))
# copy origial data
chon.no.na <- dchon
# loop over imputed variables
for(i in ch.imp.vars)
{
  print(i)
  # get the mean imputations for this variable
  ch.imp.i <- ch.imp.vars.mean[[i]]
  # get the original indices for NA
  ch.idx.i <- which(is.na(chon.no.na[[i]]))
  # replace original NA with imputed values
  chon.no.na[ch.idx.i, i] <- ch.imp.i
}




### Merge back together ############################################################
## All 3 sets
dall <- rbind(bon.no.na,sco.no.na,chon.no.na)
logs <- c("StockCode","sciname","Order","logTLinf","logTLmax","logTLBoth",
          "logWeight","logWmax","logWinf")
dimp <- dall[,logs]
dimp[,10:15] <- exp(dimp[,4:9])
names(dimp)[10:15] <- c("ImpTLinf","ImpTLmax","ImpTLBoth",
                        "ImpWeight","ImpWmax","ImpWinf")

imps <- c("StockCode","ImpTLinf","ImpTLmax","ImpTLBoth",
          "ImpWeight","ImpWmax","ImpWinf")
itraits <- dimp[,imps]
summary(itraits) 

smtraits2 <- merge(mtraits,itraits,by="StockCode",all=T)
smtraits2 <- unique(smtraits2)
write.table(smtraits2,"ICES_RAMv3.0_ALL_spp_mean_size_traits_imp_hmisc_RAMv1.0_extra.csv",sep=",",row.names=F)



