% Make a test mat file of interpolated time series from one year of COBALT
% Last year of ESM2M historical = 2005

clear all
close all

Cdir = '/Volumes/GFDL/GCM_DATA/Hindcast/';
Tdir = '/Volumes/GFDL/GCM_DATA/ESM26_hist/';

load([Cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat']);
LON  = lon;
LAT  = lat;
Z    = ncread([Cdir 'grid_spec.nc'], 'ht'); % depth

%% Water mask
[ni, nj] = size(LON);
lmask  = zeros(ni,nj);
WID = find(~isnan(temp_100(1,:,:)));
lmask(WID) = ones;

%% Retain only water cells
ID = find(Z(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = Z(ID);
GRD.lmask = lmask(ID);

%! save
save('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_ESM26_hist.mat','GRD');
        
        
        