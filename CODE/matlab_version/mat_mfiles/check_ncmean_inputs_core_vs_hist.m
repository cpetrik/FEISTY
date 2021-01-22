% check inputs, because CORE zoo biomass too low

clear all
close all

fpath = '/Volumes/MIP/GCM_DATA/CORE-forced/';
hpath = '/Volumes/FEISTY/GCM_DATA/ESM2M_hist/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Units
load([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100');
load([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb');
load([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100');
load([fpath 'ocean_cobalt_lz100_monthly_1950_2007.mat'],'lz_100');
load([fpath 'ocean_cobalt_hploss_mz100_monthly_1950_2007.mat'],'hploss_mz_100');
load([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100');
load([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat']); %,'det_btm'

%%
tp_100 = double(tp_100);
tb = double(tb);
mz_100 = double(mz_100);
lz_100 = double(lz_100);
det_btm = double(det_btm);
hploss_mz_100 = double(hploss_mz_100);
hploss_lz_100 = double(hploss_lz_100);

%% means in space
tp_mean_core = nanmean(tp_100,3);
tb_mean_core = nanmean(tb,3);
det_mean_core = nanmean(det_btm,3);
mz_mean_core = nanmean(mz_100,3);
lz_mean_core = nanmean(lz_100,3);
mzloss_mean_core = nanmean(hploss_mz_100,3);
lzloss_mean_core = nanmean(hploss_lz_100,3);

%% Compare to historical 
load([hpath 'hist_90-95_det_biom_Dmeans_Ytot.mat']); 

mz_mean_hist = double(mz_mean_hist);
lz_mean_hist = double(lz_mean_hist);
mzloss_mean_hist = double(mzloss_mean_hist);
lzloss_mean_hist = double(lzloss_mean_hist);

%%
clatlim=[-90 90];
clonlim=[-180 180];

LAT = double(geolat_t);
LON = double(geolon_t);

%% Zoo biomass
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mz_mean_core))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('CORE MZ')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mz_mean_hist))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Hist MZ')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(lz_mean_core))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('CORE LZ')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(lz_mean_hist))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('Hist LZ')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%print('-dpng',[pp 'Map_GFDL_Pre_1949_from_daily_interp_forcings.png'])

%% HP loss
figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mzloss_mean_core))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('CORE MZ loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mzloss_mean_hist))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Hist MZ loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(lzloss_mean_core))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('CORE LZ loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(lzloss_mean_hist))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('Hist LZ loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%print('-dpng',[pp 'Map_GFDL_Pre_1949_from_daily_interp_forcings.png'])

%% Det
figure(3)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_mean_core))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('CORE Det')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_mean_hist))
cmocean('tempo')
%caxis([0 2])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Hist Det')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%print('-dpng',[pp 'Map_GFDL_Pre_1949_from_daily_interp_forcings.png'])



