% Effect of diff nmort rates 
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

cmYOR=cbrewer('seq','YlOrRd',50);
close all

%% Plots in space
harv = 'All_fish03';
cfileA = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpathA=['/Volumes/GFDL/NC/Matlab_new_size/' cfileA '/'];
ppath = [pp cfileA '/'];
load([fpathA 'Means_bio_prod_fish_Climatol_All_fish03_',cfileA,'.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Afish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Aall=NaN*ones(ni,nj);
Aall(ID)=Afish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileB = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort6_BE08_noCC_RE00100';
fpathB=['/Volumes/GFDL/NC/Matlab_new_size/' cfileB '/'];
load([fpathB 'Means_bio_prod_fish_Climatol_All_fish03_',cfileB,'.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Bfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Ball=NaN*ones(ni,nj);
Ball(ID)=Bfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileC = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort7_BE08_noCC_RE00100';
fpathC=['/Volumes/GFDL/NC/Matlab_new_size/' cfileC '/'];
load([fpathC 'Means_bio_prod_fish_Climatol_All_fish03_',cfileC,'.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Cfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Call=NaN*ones(ni,nj);
Call(ID)=Cfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileD = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort4_BE08_noCC_RE00100';
fpathD=['/Volumes/GFDL/NC/Matlab_new_size/' cfileD '/'];
load([fpathD 'Means_bio_prod_fish_Climatol_All_fish03_',cfileD,'.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Dfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Dall=NaN*ones(ni,nj);
Dall(ID)=Dfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

% cfileE = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
% fpathE=['/Volumes/GFDL/NC/Matlab_new_size/' cfileE '/'];
% load([fpathE 'Means_bio_prod_fish_Climatol_All_fish03_',cfileE,'.mat'],'sf_mean','mf_mean','sp_mean',...
%     'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
% Efish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
% Eall=NaN*ones(ni,nj);
% Eall(ID)=Efish;
% clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

diffB = Ball-Aall;
diffC = Call-Aall;
diffD = Dall-Aall;

%% Subplots difference from equal red-white-blue
figure(1)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('T-dep')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffC)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('w-dep')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('T-dep, w-dep JC15')

%D
% subplot('Position',[0.5 0.1 0.5 0.4])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,diffEF)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('')
stamp(cfileA)
print('-dpng',[ppath 'Climatol_' harv '_changeAll_diff_nmort.png'])

%% LMEs
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
tlme = lme_mask_onedeg;

lme_mean = NaN*ones(66,4);
for L=1:66
    lid = find(tlme==L);
    lme_mean(L,1) = nanmean(diffB(lid));
    lme_mean(L,2) = nanmean(diffC(lid));
    lme_mean(L,3) = nanmean(diffD(lid));
%     lme_mean(L,4) = nanmean(diffE(lid));
end

lme_diffB = NaN*ones(ni,nj);
lme_diffC = lme_diffB;
lme_diffD = lme_diffB;
% lme_diffE = lme_diffB;
for L=1:66
    lid = find(tlme==L);
    lme_diffB(lid) = lme_mean(L,1);
    lme_diffC(lid) = lme_mean(L,2);
    lme_diffD(lid) = lme_mean(L,3);
%     lme_diffE(lid) = lme_mean(L,4);
end

%%
figure(2)
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('T-dep')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffC)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('w-dep')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('T-dep, w-dep JC15')

%D
% subplot('Position',[0.5 0.1 0.5 0.4])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,lme_diffEF)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('')
stamp(cfileA)
print('-dpng',[ppath 'Climatol_' harv '_changeAll_diff_nmort_LME.png'])

%%
figure(3)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Aall))
colormap(cmYOR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1.5]);
set(gcf,'renderer','painters')
title('const')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ball))
colormap(cmYOR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1.5]);
set(gcf,'renderer','painters')
title('T-dep')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Call))
colormap(cmYOR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1.5]);
set(gcf,'renderer','painters')
title('w-dep')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Dall))
colormap(cmYOR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1.5]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('T-dep, w-dep JC15')
stamp(cfileA)
print('-dpng',[ppath 'Climatol_' harv '_All_diff_nmort.png'])
