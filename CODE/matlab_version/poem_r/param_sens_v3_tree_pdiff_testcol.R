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
vec <- read.csv("Climatol_All_fish03_param_sens_v3_vecs_pdiff_lat.csv",sep=",",header = T)
pvec <- read.csv("Climatol_All_fish03_param_sens_v3_vecs_pdiff.csv",sep=",",header = T)

#raw data
vmat <- as.matrix(vec[2:39,2:6])
rvec <- as.data.frame(t(vmat))
names(rvec) <- vec$Row[2:39]

pvmat <- as.matrix(pvec[2:39,2:6])
rpvec <- as.data.frame(t(pvmat))
names(rpvec) <- pvec$Row[2:39]

## ----------------------------------- dendextend -------------------------------------------
d_iris <- dist(t(rvec)) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
dend <- as.dendrogram(hc_iris)
# Color the branches based on the clusters:
dend <- color_branches(dend, k=5) 
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
dend <- set(dend, "labels_cex", 0.5)

# And plot:
pdf(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete_pdiff_lat.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(dend, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

some_col_func <- colorspace::diverge_hcl

col_breaks = c(-15,-5,-1.1, # for low
               seq(-1,0,length=11), # for blue
               seq(0.1,1,length=10),  # for red
               1.1,5,15) # for high
test <- colorspace::diverge_hcl(26)
png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
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
png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_wardD_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
  height = 5*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_single_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_average_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_mcquitty_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_median_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_centroid_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_wardD2_heatmap_pdiff_lat.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rvec)), 
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
                  col = test,
                  breaks=col_breaks
)
dev.off()


###################    ALL PERCENT DIFFS (INCLUDING LAT)  ######################
d_fish <- dist(t(rpvec)) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "complete")
tree <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
tree <- color_branches(tree, k=5) 
# We hang the dendrogram a bit:
tree <- hang.dendrogram(tree,hang_height=0.1)
# reduce the size of the labels:
tree <- set(tree, "labels_cex", 0.5)

# And plot:
pdf(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete_pdiff.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(tree, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

some_col_func <- colorspace::diverge_hcl

col_breaks = c(-15,-5,-1.1, # for low
               seq(-1,0,length=11), # for blue
               seq(0.1,1,length=10),  # for red
               1.1,5,15) # for high
test <- colorspace::diverge_hcl(26)
png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using complete)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = tree,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()


## Diff clustering methods
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", 
                    "median", "centroid", "ward.D2")
fish_dendlist <- dendlist()
for(i in seq_along(hclust_methods)) {
  hc_fish <- hclust(d_fish, method = hclust_methods[i])   
  fish_dendlist <- dendlist(fish_dendlist, as.dendrogram(hc_fish))
}
names(fish_dendlist) <- hclust_methods



## -------------- Diff methods -------------------------------
png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_wardD_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using ward.D)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$ward.D,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_single_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using single)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$single,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_average_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using average)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$average,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_mcquitty_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using mcquitty)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$mcquitty,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_median_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using median)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$median,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_centroid_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using centroid)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$centroid,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()

png(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_wardD2_heatmap_pdiff.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
gplots::heatmap.2(as.matrix(t(rpvec)), 
                  main = "Parameter sensitivity\n (clustered using ward.D2)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = fish_dendlist$ward.D2,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Response",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  col = test,
                  breaks=col_breaks
)
dev.off()




