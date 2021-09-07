% Make GRD file for FEISTY input from 
% GFDL CORE-forced biome climatologies

clear all
close all

Cdir = '/Volumes/MIP/GCM_DATA/CORE-forced/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Depth, lat, lon, area, grid cell with seafloor
load([Cdir 'ocean_cobalt_grid.mat'])

%%  doubles
kmt = double(kmt);
depth = double(ht);
geolat_t = double(geolat_t);
geolon_t = double(geolon_t);
area = double(area_t);

%% Biomes
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);

%% Sep hemis
nid = find(geolat_t(:)>0);
sid = find(geolat_t(:)<0);
vmask = lmask;
vmask(nid) = vmask(nid)*2;
nhem = find(vmask==2);
shem = find(vmask==1);

biome8_hist = biome4_hist;
biome8_hist(nid) = biome8_hist(nid) +4;
%1 & 5: LC
%2 & 6: ECCS
%3 & 7: ECSS
%4 & 8: Coastal
pcolor(biome8_hist)
shading flat

%% Take means over each biome

for L=1:8
    lid = find(biome8_hist==L);
    Barea(L,:) = nanmean(area(lid));
    Bdepth(L,:) = nanmean(depth(lid));
    Blon(L,:) = nanmean(geolon_t(lid));
    Blat(L,:) = nanmean(geolat_t(lid));
end

%%
GRD.LON = repelem(Blon,220,1);
GRD.LAT = repelem(Blat,220,1);
GRD.Z   = repelem(Bdepth,220,1);
GRD.area = repelem(Barea,220,1);

GRD.ID = [1:length(GRD.LON)]';
GRD.N = length(GRD.LON);
GRD.lmask = ones(size(GRD.LON));

%% Save needed variables
save([Cdir 'Grid_cobalt_core_exper.mat'],'GRD');

