% Make GRD file for FEISTY input from NEMURO CCE model
% using GFDL downscaled projections

clear 
close all

Cdir = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/GFDLdown/';

%% 
load([Cdir 'feisty_gfdl_gridspec.mat'])

%Land mask
mask = BATHY;
mask(mask<10.1) = nan;
lmask = mask;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;


WID = find(lmask(:)==1);  % spatial index of water cells
NID = length(WID);

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = BATHY(ID);
GRD.lmask = lmask(ID);

%% Save needed variables
save([Cdir 'feisty_gfdl_gridspec.mat'],'mask','lmask','-append');
save([Cdir 'Data_grid_nemuro_gfdl.mat'],'GRD');
