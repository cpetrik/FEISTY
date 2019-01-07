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
library(BBmisc)

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
cfile = "Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100"
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/param_sens/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/param_sens/")

# load data
vec <- read.csv("Climatol_All_fish03_param_sens_v3_vecs.csv",sep=",",header = T)

#raw data
vmat <- as.matrix(vec[2:39,2:6])
rvec <- as.data.frame(t(vmat))
names(rvec) <- vec$Row[2:39]


#scales each column so they range -1 to 1 (uses BBmisc)
test <- normalize(vec, method = "range", range = c(-1, 1))
#I think I want to scale the whole matrix using max and min of all
test1 <- rescale(vmat, to=c(-1,1))
test1 <- as.data.frame(test1)
#same results
range01 <- function(x){ 2* (x - min(x))/(max(x)-min(x)) - 1 }
test2 <- range01(vmat)
test2 <- as.data.frame(test2)

nvec <- as.data.frame(t(test2))
names(nvec) <- vec$Row[2:39]

## ---------------------------- Ward Hierarchical Clustering -------------------------------------------
#Raw
d <- dist(t(rvec), method = "euclidean") # distance matrix
rfit <- hclust(d, method="ward.D") 
plot(rfit) # display dendogram
groups <- cutree(rfit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(rfit, k=5, border="red")

#Scaled
d <- dist(t(nvec), method = "euclidean") # distance matrix
nfit <- hclust(d)#, method="ward.D") 
plot(nfit) # display dendogram
groups <- cutree(nfit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(nfit, k=5, border="red")
#same results

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit2 <- pvclust(nvec, method.hclust="ward.D",
               method.dist="euclidean")
plot(fit2) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2, alpha=.95)

fit2 <- pvclust(nvec, method.hclust="complete",
                method.dist="euclidean")
pdf(paste0(fpath,"Climatol_All_fish03_dend_pvals.pdf"))
plot(fit2) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit2, alpha=.95)
dev.off()

## ------------------------------ Model Based Clustering -------------------------------------------
library(mclust)
fit3 <- Mclust(t(nvec))
plot(fit3) # plot results 
summary(fit3) # display the best model



## ----------------------------------- dendextend -------------------------------------------
library(dendextend)
pdf(paste0(fpath,"Climatol_All_fish03_param_sens_v3pairplot.pdf"))
pairs(t(nvec), lower.panel = NULL, cex.labels=2, pch=19, cex = 1.2)
dev.off()

d_iris <- dist(t(nvec)) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
dend <- as.dendrogram(hc_iris)

# order it the closest we can to the order of the observations:
#dend <- rotate(dend, 1:150)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=5) #, groupLabels=iris_species)

# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)

# And plot:
pdf(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(dend, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

#some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
some_col_func <- colorspace::diverge_hcl
# scaled_iris2 <- iris2 %>% as.matrix %>% scale
library(gplots)
pdf(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete_heatmap_v1.pdf"))
gplots::heatmap.2(as.matrix(t(nvec)), 
                  main = "Heatmap for the parameter sensitivity",
                  srtCol = 20,
                  dendrogram = "row",
                  Rowv = dend,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(5,0.1),      
                  key.xlab = "Cm",
                  denscol = "grey",
                  density.info = "density",
                  col = some_col_func
)
dev.off()

pdf(paste0(fpath,"Climatol_All_fish03_sensv3_hclust_complete_heatmap_v2.pdf"))
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
iris_dendlist
#correlations
iris_dendlist_cor <- cor.dendlist(iris_dendlist)
iris_dendlist_cor
library("corrplot")

pdf(paste0(fpath,"Climatol_All_fish03_sensv3_corr_diff_methods.pdf"))
corrplot::corrplot(iris_dendlist_cor, "pie", "lower")
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
#12 common sub-trees

pdf(paste0(fpath,"Climatol_All_fish03_sensv3_diff_methods.pdf"))
par(mfrow = c(4,2))
for(i in 1:8) {
  iris_dendlist[[i]] %>% set("branches_k_color", k=5) %>% plot(axes = FALSE, horiz = TRUE)
  title(names(iris_dendlist)[i])
}
dev.off()



## ----------------------------------- fn rpart -------------------------------------------
## FULL TREE 
# grow tree 
fit <- rpart(PVR ~ LongevImp + TLinf + NumSpMos + TrophImp + BatchSpawner + ReproMode + Fertilization + 
               RepGuild1 + DemersPelag + EnvTemp,
             method="anova", data=R.df)

printcp(fit) # display the results 

# visualize cross-validation results
pdf("RegressionTree/reg_tree_cross_val_PVR_rpart.pdf")
plotcp(fit) 
dev.off()

summary(fit) # detailed summary of splits
rsq.rpart(fit)	#plot approximate R-squared and relative error for different splits 
#(2 plots). labels are only appropriate for the "anova" method.

# plot tree 
pdf("RegressionTree/reg_tree_PVR_traits_rpart.pdf")
plot(fit, uniform=TRUE, 
     main="Regression Tree for PV R")
text(fit, use.n=TRUE, all=TRUE, cex=.8)
dev.off()

# create attractive postscript plot of tree 
post(fit, file = "RegressionTree/reg_tree_PVR_traits_rpart.ps", 
     title = "Regression Tree for PV R")

## PRUNED TREE 
# prune tree
#complexity parameter associated with minimum error
pfit <- prune(fit, fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
summary(pfit)

# plot the pruned tree 
pdf("RegressionTree/pruned_reg_tree_PVR_traits_rpart.pdf")
plot(pfit, uniform=TRUE, 
     main="Pruned Regression Tree for PV R")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)
dev.off()

post(pfit, file = "RegressionTree/pruned_reg_tree_PVR_traits_rpart.ps", 
     title = "Pruned Regression Tree for PV R")


