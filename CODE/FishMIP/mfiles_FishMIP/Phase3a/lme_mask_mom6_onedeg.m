% Make LME mask for MOM6 COBALT2 sims
% 1961-2010 ctrlclim & obsclim with obs fishing effort

clear all
close all

%% Map data
%cpath = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
cpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');

[ni,nj]=size(LON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%%
Pdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([Pdir 'lme_mask_om4_025_2.mat'],'lme_mask');

rpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';
load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig']);

%% diff orientation
lat = lat';
lat = fliplr(lat);
lon = lon';
lon = fliplr(lon);
% area = area';
% area = fliplr(area);
lme_mask = lme_mask';
lme_mask = fliplr(lme_mask);

%AREA_OCN = max(area,1);

%% need to shift lon
%fix lon shift
id=find(lon(:)<-180);
lon(id)=lon(id)+360;

%%
figure
pcolor(LAT); shading flat; colorbar
figure
pcolor(lat); shading flat; colorbar
figure
pcolor(LON); shading flat; colorbar
figure
pcolor(lon); shading flat; colorbar
figure
pcolor(lme_mask); shading flat; colorbar

%% regrid
tlme = griddata(lat,lon,lme_mask, LAT,LON);

figure
pcolor(tlme); shading flat; colorbar

%%
save([cpath 'lme_gfdl-mom6-cobalt2_onedeg.mat'],'tlme');%,'AREA_OCN');


