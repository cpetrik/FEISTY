% map inputs from CORE 

clear all
close all

fpath = '/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Forcing
load([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100');
load([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb');
load([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100');
load([fpath 'ocean_cobalt_lz100_monthly_1950_2007.mat'],'lz_100');
load([fpath 'ocean_cobalt_hploss_mz100_monthly_1950_2007.mat'],'hploss_mz_100');
load([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100');
load([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat']); %,'det_btm'

%% change units
tp_100 = double(tp_100);
tb = double(tb);
mz_100 = double(mz_100) * (106.0/16.0) * 12.01 * 9.0;
lz_100 = double(lz_100) * (106.0/16.0) * 12.01 * 9.0;
det_btm = double(det_btm) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
hploss_mz_100 = double(hploss_mz_100) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
hploss_lz_100 = double(hploss_lz_100) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;

%% means in space
tp_mean_core = nanmean(tp_100,3);
tb_mean_core = nanmean(tb,3);
det_mean_core = nanmean(det_btm,3);
mz_mean_core = nanmean(mz_100,3);
lz_mean_core = nanmean(lz_100,3);
mzloss_mean_core = nanmean(hploss_mz_100,3);
lzloss_mean_core = nanmean(hploss_lz_100,3);

z_mean_core = mz_mean_core + lz_mean_core;
zloss_mean_core = mzloss_mean_core + lzloss_mean_core;

%%
clatlim=[-90 90];
clonlim=[-280 80];

load coastlines

LAT = double(geolat_t);
LON = double(geolon_t);

%% 
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tp_mean_core)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_mean_core)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(z_mean_core))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_mean_core))
cmocean('tempo')
caxis([-2.5 1.5])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zloss_mean_core))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CORE_mean_forcings.png'])



