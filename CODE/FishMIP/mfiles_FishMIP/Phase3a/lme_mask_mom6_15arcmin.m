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
cdir = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([Pdir 'lme_mask_om4_025_2.mat']);
load([cdir 'esm26_area_1deg.mat']);

%% diff orientation
lat = lat';
lat = fliplr(lat);
lon = lon';
lon = fliplr(lon);
area = area';
area = fliplr(area);
lme_mask_onedeg = lme_mask_onedeg';
lme_mask_onedeg = fliplr(lme_mask_onedeg);

AREA_OCN = max(area,1);

%% need to shift lon
%fix lon shift
id=find(lon(:)>180);
lon(id)=lon(id)-360;

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
pcolor(lme_mask_onedeg); shading flat; colorbar

%% regrid
tlme = griddata(lat,lon,lme_mask_onedeg, LAT,LON);

figure
pcolor(tlme); shading flat; colorbar

%%
save([cpath 'lme_gfdl-mom6-cobalt2_onedeg_temporary.mat'],'tlme','AREA_OCN');


