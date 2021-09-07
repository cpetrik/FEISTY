% Make GRD file for FEISTY input from 
% GFDL CORE-forced

clear all
close all

Cdir = '/Volumes/MIP/GCM_DATA/CORE-forced/';

%% Depth, lat, lon, area, grid cell with seafloor
load([Cdir 'ocean_cobalt_grid.mat'])

%%  doubles
kmt = double(kmt);
depth = double(ht);
geolat_t = double(geolat_t);
geolon_t = double(geolon_t);
area = double(area_t);

%% check 
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,depth)
title('ESM2M depth')

%%
WID = find(~isnan(depth(:))); 
NID = length(WID); %48111

% Land mask
lmask = depth;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;
LID = find(lmask(:)==1);

eq1 = (WID==LID); 
sum(eq1)

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = geolon_t(ID);
GRD.LAT = geolat_t(ID);
GRD.Z   = depth(ID);
GRD.lmask = lmask(ID);
GRD.area = area(ID);

%% Save needed variables
save([Cdir 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
save([Cdir 'ocean_cobalt_grid.mat'],'lmask','-append');

