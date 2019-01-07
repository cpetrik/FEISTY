# Regression tree of parameter sensitivity test

rm(list=ls())

library(Hmisc)
#library(pastecs)
#library(fBasics)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(scatterplot3d)
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)
library(gplots)
library("corrplot")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
cfile = "Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100"
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/param_sens/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/param_sens/")

# load data
vec <- read.csv("Climatol_All_fish03_param_sens_10p_vecs_pdiff.csv",sep=",",header = T)

# magnitude
vec$mag <- sqrt(vec$F^2 + vec$P^2 + vec$D^2 + vec$Trop^2 + vec$Temp^2)
summary(vec$mag[2:39]) #1st Q = 0.018

#raw data
vmat <- as.matrix(vec[2:39,c(2:6,8)])
rvec <- as.data.frame(t(vmat))
names(rvec) <- vec$Row[2:39]


## -------------------------------- correlations ----------------------------------------
### Check that parameter changes in opposite directions have r = -1
obs1 <- names(rvec)[1:8]
obs2 <- names(rvec)[9:16]
obs3 <- names(rvec)[17:24]
obs4 <- names(rvec)[25:32]
obs5 <- names(rvec)[33:38]

pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_pairplot_response1.pdf"))
pairs((rvec[1:5,obs1]), lower.panel = panel.cor)
dev.off()
pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_pairplot_response2.pdf"))
pairs((rvec[1:5,obs2]), lower.panel = panel.cor)
dev.off()
pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_pairplot_response3.pdf"))
pairs((rvec[1:5,obs3]), lower.panel = panel.cor)
dev.off()
pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_pairplot_response4.pdf"))
pairs((rvec[1:5,obs4]), lower.panel = panel.cor)
dev.off()
pdf(paste0(fpath,"Climatol_All_fish03_param_sens10p_pairplot_response5.pdf"))
pairs((rvec[1:5,obs5]), lower.panel = panel.cor)
dev.off()

# All Corr = -1

# Get rid of params with small response (mag < 1st Q)
high <- which(rvec[6,]>=0.018)
bvec <- rvec[,c(high,20,36)]


# Multiply params with corr by -1 so both changes incr F
#h=cmax, gam=enc
newn <- c("ac-10","ac+10","ae-10","ae+10","am-10","am+10","lam-10","lam+10",
          "bc-10","bc+10","be-10","be+10","bm-10","bm+10","BE-10","BE+10",
          "fish-10","kc-10","kc+10","ke-10","ke+10","km-10","km+10","kap-10",
          "kap+10","A-10","A+10","J-10","fish+10","J+10")
names(bvec) <- newn

posF <- which(bvec[1,]>0)
negF <- which(bvec[1,]<0)
keep <- names(bvec)[posF]
mult <- names(bvec)[negF]

bvec[,30:44] <- -1*bvec[,negF]
names(bvec)[30:44] <- c("-ac+10","-ae+10","-am-10","-lam+10","-bc-10","-be-10",
                        "-bm+10","-BE-10","-fish-10","-kc+10","-ke-10","-km+10",
                        "-kap-10","-A+10","-J+10")
params <- c(posF,30:44)
tvec <- bvec[1:5,params]
#tvec <- bvec[1:5,]
tvec <- as.data.frame(t(tvec))
  
## -------------------------------- dendextend -------------------------------------------
dcorr <- as.dist((1 - cor(t(tvec)))/2)
plot(hclust(dcorr))
## Use correlations between variables "as distance"
dd <- as.dist((1 - cor(USJudgeRatings))/2)
round(1000 * dd) # (prints more nicely)
plot(hclust(dd)) # to see a dendrogram of clustered variables


d_iris <- dist(tvec) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
dend <- as.dendrogram(hc_iris)
# Color the branches based on the clusters:
dend <- color_branches(dend, k=5) 
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
dend <- set(dend, "labels_cex", 0.5)

# And plot:
pdf(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_complete_pdiff.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(dend, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

some_col_func <- colorspace::diverge_hcl

col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high
test <- colorspace::diverge_hcl(28)
png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_complete_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using complete)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = dend,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()


## Diff clustering methods
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", 
                    "median", "centroid", "ward.D2")
iris_dendlist <- dendlist()
for(i in seq_along(hclust_methods)) {
  hc_iris <- hclust(d_iris, method = hclust_methods[i])   
  iris_dendlist <- dendlist(iris_dendlist, as.dendrogram(hc_iris))
}
names(iris_dendlist) <- hclust_methods



## -------------- Diff methods -------------------------------
png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_wardD_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
  height = 5*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using ward.D)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$ward.D,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_single_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using single)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$single,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_average_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using average)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$average,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_mcquitty_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using mcquitty)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$mcquitty,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_median_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using median)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$median,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_centroid_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using centroid)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$centroid,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_wardD2_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using ward.D2)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$ward.D2,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()


##################### WARD.D LOOKS GOOD ############################
d_fish <- dist(tvec) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=5) 
# We hang the dendrogram a bit:
dfish <- hang.dendrogram(dfish,hang_height=0.1)
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.5)

# And plot:
pdf(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_wardD_pdiff.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(dfish, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high
test <- colorspace::diverge_hcl(28)
png(paste0(fpath,"Climatol_All_fish03_sens10p_hclust_wardD_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(tvec), 
                  main = "Parameter sensitivity\n (clustered using ward.D)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.8,
                  col = some_col_func #, breaks=col_breaks
)
dev.off()

## Figure out how to put strengths/magnitudes on plot
pvec <- tvec
pvec$mag <- sqrt(pvec$F^2 + pvec$P^2 + pvec$D^2 + pvec$Trop^2 + pvec$Temp^2)
write.table(pvec,"Climatol_All_fish03_param_sens_10p_table_pdiff.csv",sep=",",row.names=T)


