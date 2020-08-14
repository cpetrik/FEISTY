% Make mat files of interpolated time series from GFDL
% Preindust spinup 1850-1949

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/preindust/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-3
%tp: degC
%tb: degC

%I MAY NEED TO DIVIDE CONCENTRATIONS BY 100 m TO PUT INTO m^-2

load([fpath 'gfdl_pi_temp100_monthly_1850_1949.mat'],'temp_100');
load([fpath 'gfdl_pi_temp_btm_monthly_1850_1949.mat'],'temp_btm');
load([fpath 'gfdl_pi_zmeso100_monthly_1850_1949.mat'],'zmeso_100');
load([fpath 'gfdl_pi_det_btm_monthly_1850_1949.mat']); %,'det_btm'

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/gridspec_gfdl.mat');
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/Data_grid_gfdl.mat','GRD');

%% pcolor
geolat = double(geolat);

figure
pcolor(geolat)
shading flat
title('geolat')

figure
pcolor(geolon)
shading flat
title('geolon')

[LAT,LON] = meshgrid(lat, lon);

figure
pcolor(LAT)
shading flat
title('LAT')

figure
pcolor(LON)
shading flat
title('LON')

figure
pcolor(squeeze(temp_100(:,:,end)))
shading flat
title('T100')

%% maps
Tp=squeeze(temp_100(:,:,end));
depth = double(deptho);

clatlim=[-90 90];
clonlim=[-180 180];

elatlim=[-90 90];
elonlim=[-280 80];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,Tp)
% cmocean('thermal')
% caxis([0 35])
title('T100')

figure
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,depth)
title('GFDL Depth')

figure
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,Tp)
title('T100 on GFDL grid')

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,depth)
title('Depth on CMIP6 grid')

%% grids
% index of water cells on CMIP6 grid
[ni,nj,nt] = size(temp_100);
WID = find(~isnan(temp_100(:,:,1)));  % spatial index of water cells
NID = length(WID);                    % number of water cells

% index of water cells of GFDL grid
[mi,mj] = size(deptho);
GID = find(~isnan(deptho(:)));        % spatial index of water cells
MID = length(GID);                    % number of water cells

% WID NOT EQUAL TO GRD.ID!!!

% indexes
[m,n] = ind2sub([ni,nj],WID(1));
[j,k] = ind2sub([mi,mj],MID(1));

tp = temp_100(m,n,1);
tb = temp_btm(m,n,1);
d = det_btm(m,n,1);
z = zmeso_100(m,n,1);
ilon = LON(m,n);
ilat = LAT(m,n);

hc = depth(m,n,1);
clon = geolon(m,n);
clat = geolat(m,n);
he = depth(j,k);
elon = geolon(j,k);
elat = geolat(j,k);

% DAILY INTERPOLATIONS CORRECT, GRID INFO NOT MATCHED UP WITH IT

%% regrid depth to cmip6 lat and lon
clear depth
deptho = double(deptho);

% geolon(:,1:12) = repmat(geolon(:,13),1,12);
% geolat(:,1:12) = repmat([-89.5:-78.5],360,1);
% deptho(isnan(deptho(:))) = 0;
% 
% depth = griddata(deptho,geolat,geolon,LAT,LON); %DID NOT WORK!

% figure
% axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat,geolon,deptho)
% title('GFDL Depth')
% 
% figure
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,depth)
% title('CMIP6 depth')

%%
mlon=find(geolon(:,end)<-179);
depth(1:180,:) = deptho(mlon:end,:);
depth(181:360,:) = deptho(1:(mlon-1),:);

figure
axesm ('Robinson','MapLatLimit',elatlim,'MapLonLimit',elonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,deptho)
title('GFDL Depth')

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,depth)
title('CMIP6 depth')

DID = find(~isnan(depth(:))); 
IDN = length(DID);
h = depth(m,n,1);
