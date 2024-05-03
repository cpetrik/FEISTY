% Make GRD file for FEISTY input from NEMURO CCE model
% using IPSL downscaled projections

clear 
close all

Cdir = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
%Cdir = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/ESM_data/NEMURO/';
Sdir = '/Users/cpetrik/Documents/NEMURO/';

%% 
load([Cdir 'feisty_ipsl_gridspec.mat'])

%Land mask
mask = BATHY;
mask(mask<10.1) = nan;
lmask = mask;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;


WID = find(lmask(:)==1);  % spatial index of water cells
NID = length(WID);

%% Area calc
%It is a regular 0.1x0.1 deg grid, so you can easily calculate the area (in km2) of each grid cell as:
%area = dx * dy = (0.1*111.1*cos(lat)) * (0.1*111.1) = 123.43*cos(lat)
%Area = abs(123.43*cos(LAT));

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = BATHY(ID);
GRD.lmask = lmask(ID);

%% Save needed variables
save([Cdir 'feisty_ipsl_gridspec.mat'],'mask','lmask','-append');
save([Cdir 'Data_grid_nemuro_ipsl.mat'],'GRD');
%save([Sdir 'Data_grid_nemuro_ipsl.mat'],'GRD');
