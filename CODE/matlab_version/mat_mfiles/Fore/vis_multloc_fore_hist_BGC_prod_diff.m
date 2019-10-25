% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100
% production of phyto, zoop, and fish

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%colors
cmPO=cbrewer('div','PuOr',50,'PCHIP');
cmPO = flipud(cmPO);

%% NPP and zoop
load([bpath 'cobalt_det_temp_zoop_npp_means.mat']);

% molN/m2/s --> g/m2/d
mzloss_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_hist = ptemp_mean_hist - 273;

mzloss_fore = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_fore = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_fore = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
ptemp_fore = ptemp_mean_fore - 273;

zloss_hist = mzloss_hist + lzloss_hist;
zloss_fore = mzloss_fore + lzloss_fore;

ZlDet_hist = zloss_hist./det_hist;
ZlDet_fore = zloss_fore./det_fore;
lZlDet_hist = log10(zloss_hist./det_hist);
lZlDet_fore = log10(zloss_fore./det_fore);

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%% Percent differences
diffN = (npp_fore-npp_hist) ./ npp_hist;
diffMZ = (mzloss_fore-mzloss_hist) ./ mzloss_hist;
diffLZ = (lzloss_fore-lzloss_hist) ./ lzloss_hist;
diffZ = (zloss_fore-zloss_hist) ./ zloss_hist;
diffDet = (det_fore-det_hist) ./ det_hist;
diffZD = (ZlDet_fore-ZlDet_hist) ./ ZlDet_hist;
diffLZD = (lZlDet_fore-lZlDet_hist) ./ lZlDet_hist;
diffPT = (ptemp_fore-ptemp_hist);

%% Maps
% All 4 on subplots
figure(1)
% Temp
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffPT))
colormap(cmPO)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-7 7]);
colorbar('Position',[0.1 0.525 0.3 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
%title('Forecast - Hindcast Temp')
title('Temperature')

% NPP
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffN)*100)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-75 75]);
set(gcf,'renderer','painters')
%title('Forecast - Hindcast Detritus')
title('NPP')

% Det
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffDet)*100)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-75 75]);
colorbar('Position',[0.6 0.525 0.3 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
%title('Forecast - Hindcast NPP')
title('Bottom Detrital Flux')

% Zoop
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(diffZ)*100)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-75 75]);
set(gcf,'renderer','painters')
%title('Forecast - Hindcast Mesozoo')
title('Mesozoo')
%stamp(cfile)
print('-dpng',[pp 'Hist_Fore_global_prod_diff_temp_BGC.png'])

