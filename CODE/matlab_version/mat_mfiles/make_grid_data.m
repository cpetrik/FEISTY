% Make a test mat file of interpolated time series from one year of COBALT
% Last year of ESM2M historical = 2005

clear all
close all

Cdir = '/Volumes/GFDL/GCM_DATA/Hindcast/';
Tdir = '/Volumes/GFDL/GCM_DATA/ESM2M_hist/';

TIME = ncread([Tdir 'ocean_cobalt_biomass_100.186101-200512.nlgz_100.nc'],'average_T1'); 
LON  = ncread([Cdir 'grid_spec.nc'], 'geolon_t'); % lon
LAT  = ncread([Cdir 'grid_spec.nc'], 'geolat_t'); % lon
Z    = ncread([Cdir 'grid_spec.nc'], 'ht'); % depth
zt   = ncread([Cdir 'grid_spec.nc'], 'zt'); % depth levels
dx   = ncread([Cdir 'grid_spec.nc'], 'dxt');
dy   = ncread([Cdir 'grid_spec.nc'], 'dyt');
kmt  = ncread([Cdir 'grid_spec.nc'], 'kmt');
dxtn = ncread([Cdir 'grid_spec.nc'], 'dxtn');
dyte = ncread([Cdir 'grid_spec.nc'], 'dyte');
dat  = dx.*dy; % area in m
datr = 1.0./(dat+eps);

%% Water mask
[ni, nj] = size(LON);
nk     = length(zt);
lmask  = zeros(ni,nj,nk);
 
for k=1:nk
	 for j=1:nj
			for i=1:ni
				 if (kmt(i,j) >= k)
						lmask(i,j,k) = 1.0;
				 else
						lmask(i,j,k) = 0.0;
				 end
			end
	 end
end

lmask = lmask(:,:,1);

%% Retain only water cells
ID = find(Z(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = Z(ID);
GRD.DX = dx(ID);
GRD.DY = dy(ID);
GRD.AREA  = dat(ID);
GRD.dxtn  = dxtn(ID);
GRD.dyte  = dyte(ID);
GRD.datr  = datr(ID);
GRD.lmask = lmask(ID);

%! save
save('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat','GRD');
        
        
        