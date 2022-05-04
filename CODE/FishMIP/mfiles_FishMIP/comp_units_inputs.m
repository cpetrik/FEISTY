% Compare COBALT values and ranges
% Climatol (mg C), Hist/CORE (mol N), for 0.25 deg (mol C)

clear all
close all

%% Agreggated
% load('cobalt_det_biom_means.mat','geolon_t','geolat_t','lon','lat',...
%     'det_mean_hist','det_mean_clim','det_5yr_hist','det_1yr_hist');
% 
% load('cobalt_npp_means.mat',...
%     'npp_mean_hist','npp_mean_clim','npp_5yr_hist','npp_1yr_hist','mo_hist');
% 
% load('cobalt_zoop_biom_means.mat',...
%     'mz_mean_hist','lz_mean_hist','mzloss_mean_hist','lzloss_mean_hist',...
%     'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim',...
%     'mz_5yr_hist','lz_5yr_hist');

load('cobalt_det_temp_zoop_npp_means.mat')%,'det_mean_hist','npp_mean_hist','mz_mean_hist','lz_mean_hist');

%% Climatol
% cpath='/Volumes/MIP/GCM_DATA/ESM26_hist/';
% 
% % in orig units of mgC/m2 or mgC/m2/d
% load([cpath 'npp_100_1deg_ESM26_5yr_clim_191_195.mat'])
% load([cpath 'fcdet_btm_1deg_ESM26_5yr_clim_191_195.mat'])
% load([cpath 'mdz_100_1deg_ESM26_5yr_clim_191_195.mat'])
% load([cpath 'lgz_100_1deg_ESM26_5yr_clim_191_195.mat'])
% 
% npp_100(npp_100<0) = 0;
% mdz_100(mdz_100<0) = 0;
% lgz_100(lgz_100<0) = 0;
% 
%% Hist
% hpath='/Volumes/MIP/GCM_DATA/ESM2M_hist/';
% 
% load('cobalt_hist_npp_zprod_det_nanmeans.mat','det_mean_hist','npp_mean_hist',...
%     'det_ann_hist','npp_ann_hist','det_5yr_hist','npp_5yr_hist',...
%     'det_1yr_hist','npp_1yr_hist');


%% Conversions
%molN/m2/s --> molC/m2/s
%106/16 mol C in 1 mol N
det_hist = det_mean_hist * (106.0/16.0);
npp_hist = npp_mean_hist * (106.0/16.0);

%molN/m2/s --> molC/m3
%106/16 mol C in 1 mol N
mz_hist = mz_mean_hist * (106.0/16.0) * (1/100);
lz_hist = lz_mean_hist * (106.0/16.0) * (1/100);

%molN/m2/s --> molC/m2
%106/16 mol C in 1 mol N
%100 m integrated
mz_vint_hist = mz_mean_hist * (106.0/16.0);
lz_vint_hist = lz_mean_hist * (106.0/16.0);


%mgC/m2/d --> molC/m2/s
%12.01 g C in 1 mol C
%1e3 mg in 1 g 
det_clim = det_mean_clim * 1e-3 * (1/12.01) * 1/(60*60*24);
npp_clim = npp_mean_clim * 1e-3 * (1/12.01) * 1/(60*60*24);

%mgC/m2/s --> molC/m3
%12.01 g C in 1 mol C
%1e3 mg in 1 g
%100 m integrated
mz_clim = mz_mean_clim * 1e-3 * (1/12.01) * (1/100);
lz_clim = lz_mean_clim * 1e-3 * (1/12.01) * (1/100);

%mgC/m2/s --> molC/m2
%12.01 g C in 1 mol C
%1e3 mg in 1 g
%100 m integrated
mz_vint_clim = mz_mean_clim * 1e-3 * (1/12.01);
lz_vint_clim = lz_mean_clim * 1e-3 * (1/12.01);

%% ranges
Cmaxmin(1,1) = nanmax(npp_clim(:));
Cmaxmin(2,1) = nanmax(det_clim(:));
Cmaxmin(3,1) = nanmax(mz_clim(:)+lz_clim(:));
Cmaxmin(4,1) = nanmax(mz_vint_clim(:)+lz_vint_clim(:));
Cmaxmin(1,2) = nanmin(npp_clim(:));
Cmaxmin(2,2) = nanmin(det_clim(:));
Cmaxmin(3,2) = nanmin(mz_clim(:)+lz_clim(:));
Cmaxmin(4,2) = nanmin(mz_vint_clim(:)+lz_vint_clim(:));

Hmaxmin(1,1) = nanmax(npp_hist(:));
Hmaxmin(2,1) = nanmax(det_hist(:));
Hmaxmin(3,1) = nanmax(mz_hist(:)+lz_hist(:));
Hmaxmin(4,1) = nanmax(mz_vint_hist(:)+lz_vint_hist(:));
Hmaxmin(1,2) = nanmin(npp_hist(:));
Hmaxmin(2,2) = nanmin(det_hist(:));
Hmaxmin(3,2) = nanmin(mz_hist(:)+lz_hist(:));
Hmaxmin(4,2) = nanmin(mz_vint_hist(:)+lz_vint_hist(:));

%0.25 deg
maxmin(1,1) = 2.7e-6;
maxmin(2,1) = 2.6e-6;
maxmin(3,1) = 0.008;
maxmin(4,1) = 0.3;
maxmin(1,2) = -1.95e-15; %Need to remove negatives?
maxmin(2,2) = 5e-17;
maxmin(3,2) = 1.5e-9;
maxmin(4,2) = 3.3e-8;

%% maps
clatlim=[-90 90];
clonlim=[-180 180];

geolat_t = double(geolat_t);
geolon_t = double(geolon_t);

%% Clim
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat,lon,npp_clim)
cmocean('tempo')
caxis([0 2e-6])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('Clim NPP-vint')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat,lon,det_clim)
cmocean('tempo')
caxis([5e-14 1.4e-7])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Clim Expc-bot')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat,lon,mz_clim+lz_clim)
cmocean('tempo')
caxis([1.5e-9 0.004])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('Clim Zmeso')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat,lon,mz_vint_clim+lz_vint_clim)
cmocean('tempo')
caxis([3e-8 0.1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('Clim Zmeso-vint')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',['Map_Clim_det_npp_zoop_forcing_means_molC.png'])


%% Hist
figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,npp_hist)
cmocean('tempo')
caxis([0 2e-6])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('Hist NPP')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,det_hist)
cmocean('tempo')
caxis([5e-14 1.4e-7])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Hist Expc-bot')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,mz_hist+lz_hist)
cmocean('tempo')
caxis([1.5e-9 0.004])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('Hist Zmeso')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,mz_vint_hist+lz_vint_hist)
cmocean('tempo')
caxis([3e-8 0.1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('Hist Zmeso-vint')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',['Map_Hist_det_npp_zoop_forcing_means_molC.png'])

%% chl, MLD from CMIP6
% phyc, zooc, thetao, thkcello,
mpath = '/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist';

load('gfdl_hist_surf_chl_monthly_1950_2014.mat');
load('gfdl_hist_mld_monthly_1965_2014.mat');
load('gfdl_hist_npp_monthly_1951_2014.mat');
load('.mat');











