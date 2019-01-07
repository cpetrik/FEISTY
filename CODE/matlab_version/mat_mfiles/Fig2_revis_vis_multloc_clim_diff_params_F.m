% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/v4_kt85_BE75/';
cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_area_1deg.mat']);
AREA_OCN = max(area,1);

%
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

close all

% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


% colors
cmYGB=cbrewer('seq','YlGnBu',50);
cmPBG=cbrewer('seq','PuBuGn',50);
cmPB=cbrewer('seq','BuPu',50);

%% Plots in space

cfileA = 'Dc_enc70-b200_m4-b250-k086_c20-b250_D075_J100_A100_Sm025_nmort1_BE08_noCC_RE00100';
fpathA=['/Volumes/GFDL/NC/Matlab_new_size/' cfileA '/'];
load([fpathA 'Means_bio_prod_fish_Climatol_' harv '_' cfileA '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Asf=NaN*ones(ni,nj);
Amf=NaN*ones(ni,nj);
Asf(ID)=sf_mean;
Amf(ID)=mf_mean;
Asp=NaN*ones(ni,nj);
Amp=NaN*ones(ni,nj);
Alp=NaN*ones(ni,nj);
Asp(ID)=sp_mean;
Amp(ID)=mp_mean;
Alp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

cfileB = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A100_Sm025_nmort1_BE08_noCC_RE00100';
fpathB=['/Volumes/GFDL/NC/Matlab_new_size/' cfileB '/'];
load([fpathB 'Means_bio_prod_fish_Climatol_' harv '_' cfileB '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Bsf=NaN*ones(ni,nj);
Bmf=NaN*ones(ni,nj);
Bsf(ID)=sf_mean;
Bmf(ID)=mf_mean;
Bsp=NaN*ones(ni,nj);
Bmp=NaN*ones(ni,nj);
Blp=NaN*ones(ni,nj);
Bsp(ID)=sp_mean;
Bmp(ID)=mp_mean;
Blp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

cfileC = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathC=['/Volumes/GFDL/NC/Matlab_new_size/' cfileC '/'];
ppath = [pp cfileC '/'];
load([fpathC 'Means_bio_prod_fish_Climatol_' harv '_' cfileC '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Csf=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Csf(ID)=sf_mean;
Cmf(ID)=mf_mean;
Csp=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Csp(ID)=sp_mean;
Cmp(ID)=mp_mean;
Clp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

% Diff maps of all fish
AF = Asf+Amf;
BF = Bsf+Bmf;
CF = Csf+Cmf;
AP = Asp+Amp+Alp;
BP = Bsp+Bmp+Blp;
CP = Csp+Cmp+Clp;

%% 
figure(1)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AF))
colormap(cmYGB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-1.0,1.55,'Forage Fishes')

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AP))
colormap(cmYGB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.0,1.55,'Large Pelagics')

%B
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BF))
colormap(cmYGB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BP))
colormap(cmYGB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')

%C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CF))
colormap(cmYGB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
%stamp([harv '_' cfile])

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CP))
colormap(cmYGB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.04 0.34 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
print('-dpng',[ppath 'Fig2_Climatol_' harv '_FPcomp_params_YGB.png'])

%%
figure(2)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AF))
colormap(cmPBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-1.0,1.55,'Forage Fishes')

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AP))
colormap(cmPBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.0,1.55,'Large Pelagics')

%B
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BF))
colormap(cmPBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BP))
colormap(cmPBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')

%C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CF))
colormap(cmPBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
%stamp([harv '_' cfile])

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CP))
colormap(cmPBG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.04 0.34 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
print('-dpng',[ppath 'Fig2_Climatol_' harv '_FPcomp_params_PBG.png'])

%%
figure(3)
%A
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AF))
colormap(cmPB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')
text(-1.0,1.55,'Forage Fishes')

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AP))
colormap(cmPB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.0,1.55,'Large Pelagics')

%B
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BF))
colormap(cmPB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(BP))
colormap(cmPB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')

%C
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CF))
colormap(cmPB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')
%stamp([harv '_' cfile])

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(CP))
colormap(cmPB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.04 0.34 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
print('-dpng',[ppath 'Fig2_Climatol_' harv '_FPcomp_params_PB.png'])

