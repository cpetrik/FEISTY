% Map depth at 3 locations
% ESM2.6 Climatology of 5 yrs

clear 
close all

%% grid
cpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
Pdir = '/Volumes/petrik-lab/Feisty/GCM_Data/ESM26_hist/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

load coastlines;                     %decent looking coastlines

%% Ind plots
figure(1)
%  All Pac
axesm ('eqdcylin','MapLatLimit',[-30 80],'MapLonLimit',[-200 -60],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,depth)
cmocean('deep')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2 3]);
colorbar
set(gcf,'renderer','painters')
%title('log10 mean NPP (mgC m^-^2)')
print('-dpng',[pp 'Climatol_map_depth_3locs.png'])

%% EBS
figure(2)
axesm ('eqdcylin','MapLatLimit',[45 75],'MapLonLimit',[-200 -145],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,depth)
cmocean('deep')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2 3]);
colorbar
set(gcf,'renderer','painters')
%title('log10 mean NPP (mgC m^-^2)')
print('-dpng',[pp 'Climatol_map_depth_EBS.png'])

%% HOT
figure(5)
axesm ('eqdcylin','MapLatLimit',[0 45],'MapLonLimit',[-200 -110],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,depth)
cmocean('deep')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2 3]);
colorbar
set(gcf,'renderer','painters')
%title('log10 mean NPP (mgC m^-^2)')
print('-dpng',[pp 'Climatol_map_depth_HOT.png'])

%% Peru Up
figure(6)
axesm ('eqdcylin','MapLatLimit',[-30 20],'MapLonLimit',[-150 -60],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,depth)
cmocean('deep')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([2 3]);
colorbar
set(gcf,'renderer','painters')
%title('log10 mean NPP (mgC m^-^2)')
print('-dpng',[pp 'Climatol_map_depth_PUP.png'])


