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
vec <- read.csv("Climatol_All_fish03_param_sens_v2_vecs_pdiff.csv",sep=",",header = T)

#raw data
vmat <- as.matrix(vec[2:39,2:6])
rvec <- as.data.frame(t(vmat))
names(rvec) <- vec$Row[2:39]

#scale so all ranges -1 to 1 
range01 <- function(x){ 2* (x - min(x))/(max(x)-min(x)) - 1 }
test2 <- range01(vmat)
test2 <- as.data.frame(test2)
#all become -1, don't use
nvec <- as.data.frame(t(test2))
names(nvec) <- vec$Row[2:39]

## ---------------------------- Ward Hierarchical Clustering -------------------------------------------
#Ward
d <- dist(t(rvec), method = "euclidean") # distance matrix
rfit <- hclust(d, method="ward.D") 
plot(rfit) # display dendogram
groups <- cutree(rfit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(rfit, k=5, border="red")

#Complete
d <- dist(t(rvec), method = "euclidean") # distance matrix
rfit <- hclust(d, method="complete") 
plot(rfit) # display dendogram
groups <- cutree(rfit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(rfit, k=5, border="red")


# Ward Hierarchical Clustering with Bootstrapped p values
fit2 <- pvclust(rvec, method.hclust="ward.D",
               method.dist="euclidean")
pdf(paste0(fpath,"Climatol_All_fish03_pvclust_wardD_pvals_pdiff.pdf"))
plot(fit2) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2, alpha=.95)
dev.off()

fit3 <- pvclust(rvec, method.hclust="complete",
                method.dist="euclidean")
pdf(paste0(fpath,"Climatol_All_fish03_pvclust_complete_pvals_pdiff.pdf"))
plot(fit3) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit3, alpha=.95)
dev.off()



## ----------------------------------- dendextend -------------------------------------------
pdf(paste0(fpath,"Climatol_All_fish03_param_sens_pairplot_pdiff.pdf"))
pairs(t(nvec), lower.panel = NULL, cex.labels=2, pch=19, cex = 1.2)
dev.off()

d_iris <- dist(t(nvec)) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
dend <- as.dendrogram(hc_iris)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=5) #, groupLabels=iris_species)

# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)

# And plot:
pdf(paste0(fpath,"Climatol_All_fish03_hclust_complete_pdiff.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(dend, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

#some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
some_col_func <- colorspace::diverge_hcl
# scaled_iris2 <- iris2 %>% as.matrix %>% scale
pdf(paste0(fpath,"Climatol_All_fish03_hclust_complete_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using complete)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = dend,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
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
#correlations
iris_dendlist_cor <- cor.dendlist(iris_dendlist)
iris_dendlist_cor

pdf(paste0(fpath,"Climatol_All_fish03_corr_diff_methods_pdiff.pdf"))
corrplot(iris_dendlist_cor, "pie", "lower")
dev.off()

#compare
# The `which` parameter allows us to pick the elements in the list to compare
iris_dendlist %>% dendlist(which = c(1,3)) %>% ladderize %>% 
  set("branches_k_color", k=3) %>% 
  # untangle(method = "step1side", k_seq = 3:20) %>%
  tanglegram(faster = TRUE) # (common_subtrees_color_branches = TRUE)

# The `which` parameter allows us to pick the elements in the list to compare
iris_dendlist %>% dendlist(which = c(1,3)) %>% ladderize %>% 
  # untangle(method = "step1side", k_seq = 3:20) %>%
  set("rank_branches") %>%
  tanglegram(common_subtrees_color_branches = TRUE)
length(unique(common_subtrees_clusters(iris_dendlist[[1]], iris_dendlist[[3]]))[-1])
#9 common sub-trees

pdf(paste0(fpath,"Climatol_All_fish03_diff_methods_pdiff.pdf"))
par(mfrow = c(4,2))
for(i in 1:8) {
  iris_dendlist[[i]] %>% set("branches_k_color", k=5) %>% plot(axes = FALSE, horiz = TRUE)
  title(names(iris_dendlist)[i])
}
dev.off()



## -------------- Diff methods -------------------------------
pdf(paste0(fpath,"Climatol_All_fish03_hclust_wardD_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using ward.D)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$ward.D,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_hclust_single_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using single)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$single,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_hclust_average_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using average)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$average,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_hclust_mcquitty_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using mcquitty)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$mcquitty,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_hclust_median_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using median)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$median,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_hclust_centroid_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using centroid)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$centroid,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_hclust_wardD2_heatmap_pdiff.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Parameter sensitivity\n (clustered using ward.D2)",
                  srtCol = 60,
                  dendrogram = "row",
                  Rowv = iris_dendlist$ward.D2,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Normalized response",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()
