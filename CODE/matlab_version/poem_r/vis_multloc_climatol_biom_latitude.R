################################################################################

# Boxplots of biomass vs. latitude

################################################################################

rm(list=ls())

library(lattice) #multipanel plots
library(gridExtra) #gridarrange
library(ggplot2) 

#PU laptop
source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
cfile = "Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100"
setwd(paste0("/Volumes/GFDL/NC/Matlab_new_size/", cfile, "/"))
fpath <- paste0("/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/", cfile, "/")

# load data
Lat <- read.csv(paste0("Lat_bio_Climatol_All_fish03_",cfile,".csv"),sep=",",header = T,stringsAsFactors = T)
Lat$lat_bin <- as.factor(Lat$lat_bin)

Lat$depth <- "Basin"
Lat$depth[Lat$coast==1] <- "Shelf"
Lat$depth <- as.factor(Lat$depth)
Lat$depth <- relevel(Lat$depth, "Shelf")

Lat$AllF10 <- log10(Lat$AllF)
Lat$AllP10 <- log10(Lat$AllP)
Lat$AllD10 <- log10(Lat$AllD)
Lat$All10 <- log10(Lat$All)
Lat$AllS10 <- log10(Lat$AllS)
Lat$AllM10 <- log10(Lat$AllM)
Lat$AllL10 <- log10(Lat$AllL)

#scale_fill_grey()
#scale_fill_brewer(palette="Paired")
#stat_summary(fun.y=median, geom="line", aes(group=1), colour="blue")
  
#Biomass Box plots
xF<-ggplot(Lat, aes(x = lat_bin, y = AllF10, fill=depth)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) +
  theme(text = element_text(size=12), legend.position = c(0.2,0.5), 
        legend.title = element_blank()) + 
  ylab("") + scale_fill_brewer(palette="Paired") + 
  scale_colour_brewer(palette="Paired") +
  stat_summary(fun.y=median, geom="line", aes(group=depth)) + 
  stat_summary(fun.y="median", geom="point", colour="black", size=1) +
  aes(colour = depth) +
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="A. Forage fish") + 
  scale_y_continuous(breaks = seq(-4, 2, 1),limits=c(-4, 2)) + 
  scale_x_discrete(name = "Latitude",
                   labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))

xP<-ggplot(Lat, aes(x = lat_bin, y = AllP10, fill=depth)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) +
  theme(text = element_text(size=12), legend.position = "none") + 
  ylab("") + scale_fill_brewer(palette="Paired") + 
  scale_colour_brewer(palette="Paired") +
  stat_summary(fun.y=median, geom="line", aes(group=depth)) + 
  stat_summary(fun.y="median", geom="point", colour="black", size=1) +
  aes(colour = depth) +
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="B. Large pelagics") + 
  scale_y_continuous(breaks = seq(-4, 2, 1),limits=c(-4, 2)) + 
  scale_x_discrete(name = "",labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))

xD<-ggplot(Lat, aes(x = lat_bin, y = AllD10, fill=depth)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) +
  theme(text = element_text(size=12), legend.position = "none") + 
  ylab("") + scale_fill_brewer(palette="Paired") + 
  scale_colour_brewer(palette="Paired") +
  stat_summary(fun.y=median, geom="line", aes(group=depth)) + 
  stat_summary(fun.y="median", geom="point", colour="black", size=1) +
  aes(colour = depth) +
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="C. Demersals",y = expression(biomass~(g~m^-2))) + 
  scale_y_continuous(breaks = seq(-4, 2, 1),limits=c(-4, 2)) + 
  scale_x_discrete(name = "Latitude",
                   labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))

xA<-ggplot(Lat, aes(x = lat_bin, y = All10, fill=depth)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) +
  theme(text = element_text(size=12), legend.position = "none") + 
  ylab("") + scale_fill_brewer(palette="Paired") + 
  scale_colour_brewer(palette="Paired") +
  stat_summary(fun.y=median, geom="line", aes(group=depth)) + 
  stat_summary(fun.y="median", geom="point", colour="black", size=1) +
  aes(colour = depth) +
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="D. All fishes",y = expression(biomass~(g~m^-2))) + 
  scale_y_continuous(breaks = seq(-4, 2, 1),limits=c(-4, 2)) + 
  scale_x_discrete(name = "",labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))

pdf(paste0(fpath,"Climatol_All_fish03_biom_latitude_depth_boxplot.pdf"))
grid.arrange(xF,xP,xD,xA,ncol=2)
dev.off()

# Fractions
xPD<-ggplot(Lat, aes(x = lat_bin, y = FracPD)) + theme_bw() + geom_boxplot(fill='grey80',outlier.size=0.25) +
  theme(text = element_text(size=12)) + ylab("") + 
  stat_summary(fun.y=median, geom="line", aes(group=1), colour="blue") + 
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="A", y = expression(P/(P+D))) + 
  scale_x_discrete(name = "Latitude",
                   labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))
xPF<-ggplot(Lat, aes(x = lat_bin, y = FracPF)) + theme_bw() + geom_boxplot(fill='grey80',outlier.size=0.25) +
  theme(text = element_text(size=12)) + ylab("") + 
  stat_summary(fun.y=median, geom="line", aes(group=1), colour="blue") + 
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="B", y = expression(P/(P+F))) + 
  scale_x_discrete(labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))
xLM<-ggplot(Lat, aes(x = lat_bin, y = FracLM)) + theme_bw() + geom_boxplot(fill='grey80',outlier.size=0.25) +
  theme(text = element_text(size=12)) + ylab("") + 
  stat_summary(fun.y=median, geom="line", aes(group=1), colour="blue") + 
  coord_flip() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title="C", y = expression(L/(L+M))) + 
  scale_x_discrete(name = "Latitude",
                   labels=c("", "-80", "","-70","","-60","", "-50", "","-40","","-30","", "-20", "","-10","","0","", "10", "","20","","30","", "40", "","50","","60","", "70", "","80",""))

pdf(paste0(fpath,"Climatol_All_fish03_fracs_latitude_boxplot.pdf"))
grid.arrange(xPD,xPF,xLM,ncol=2)
dev.off()



