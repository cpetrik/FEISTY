% Effect of diff harvest rates on F vs. Large
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

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
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in space
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

harv = 'All_fish03';
load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat'],...
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

harvB = 'fish_F030_P060_D060';
load([fpath 'Means_bio_prod_fish_Climatol_' harvB '_' cfile '.mat'],...
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

harvC = 'fish_F030_P015_D015';
load([fpath 'Means_bio_prod_fish_Climatol_' harvC '_' cfile '.mat'],...
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

harvD = 'fish_F015_P030_D030';
load([fpath 'Means_bio_prod_fish_Climatol_' harvD '_' cfile '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Dsf=NaN*ones(ni,nj);
Dmf=NaN*ones(ni,nj);
Dsf(ID)=sf_mean;
Dmf(ID)=mf_mean;
Dsp=NaN*ones(ni,nj);
Dmp=NaN*ones(ni,nj);
Dlp=NaN*ones(ni,nj);
Dsp(ID)=sp_mean;
Dmp(ID)=mp_mean;
Dlp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

harvE = 'fish_F060_P030_D030';
load([fpath 'Means_bio_prod_fish_Climatol_' harvE '_' cfile '.mat'],...
    'sf_mean','mf_mean','sp_mean','mp_mean','lp_mean');
Esf=NaN*ones(ni,nj);
Emf=NaN*ones(ni,nj);
Esf(ID)=sf_mean;
Emf(ID)=mf_mean;
Esp=NaN*ones(ni,nj);
Emp=NaN*ones(ni,nj);
Elp=NaN*ones(ni,nj);
Esp(ID)=sp_mean;
Emp(ID)=mp_mean;
Elp(ID)=lp_mean;
clear sf_mean mf_mean sp_mean mp_mean lp_mean

% Diff maps of all fish
AF = Asf+Amf;
BF = Bsf+Bmf;
CF = Csf+Cmf;
DF = Dsf+Dmf;
EF = Esf+Emf;
AP = Asp+Amp+Alp;
BP = Bsp+Bmp+Blp;
CP = Csp+Cmp+Clp;
DP = Dsp+Dmp+Dlp;
EP = Esp+Emp+Elp;

diffBF = BF-AF;
diffCF = CF-AF;
diffDF = DF-AF;
diffEF = EF-AF;

%% Subplots difference from equal red-white-blue
figure(1)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffBF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('F030, P060, D060')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffCF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('F030, P015, D015')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffDF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('F015, P030, D030')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffEF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F060, P030, D030')
stamp(cfile)
print('-dpng',[ppath 'Climatol_' harv '_changeF_diffHarv.png'])

%% LMEs
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/'; 
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
tlme = lme_mask_onedeg;

lme_mean = NaN*ones(66,4);
for L=1:66
    lid = find(tlme==L);
    lme_mean(L,1) = nanmean(diffBF(lid));
    lme_mean(L,2) = nanmean(diffCF(lid));
    lme_mean(L,3) = nanmean(diffDF(lid));
    lme_mean(L,4) = nanmean(diffEF(lid));
end

lme_diffBF = NaN*ones(ni,nj);
lme_diffCF = lme_diffBF;
lme_diffDF = lme_diffBF;
lme_diffEF = lme_diffBF;
for L=1:66
    lid = find(tlme==L);
    lme_diffBF(lid) = lme_mean(L,1);
    lme_diffCF(lid) = lme_mean(L,2);
    lme_diffDF(lid) = lme_mean(L,3);
    lme_diffEF(lid) = lme_mean(L,4);
end

%%
figure(2)
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffBF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('F030, P060, D060')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffCF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('F030, P015, D015')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffDF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('F015, P030, D030')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffEF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('F060, P030, D030')
stamp(cfile)
print('-dpng',[ppath 'Climatol_' harv '_changeF_diffHarv_LME.png'])
