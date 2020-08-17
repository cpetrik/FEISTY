% Make GRD file for FEISTY input from CESM 1 degree model

clear all
close all

Cdir = '/Volumes/FEISTY/Fish-MIP/CMIP6/';

%% Lat & Lon
ncid = netcdf.open([Cdir 'IPSL/hist/ipsl-cm6a-lr_r1i1p1f1_historical_tob_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
%%
tob(tob >= 1.00e+20) = NaN;

%Land mask
mask = squeeze(tob(:,:,1));
lmask = double(mask);
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);

LID = find(lmask(:)==1);
WID = find(~isnan(tob(:,:,1)));  % spatial index of water cells
NID = length(WID);

eq1 = (WID==LID); %41328
sum(eq1)

%% Depth, lat, lon
ncdisp([Cdir 'IPSL/mesh_mask_eORCA1_v2.2.nc'])

ncid = netcdf.open([Cdir 'IPSL/mesh_mask_eORCA1_v2.2.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%% seafloor depths
%mbathy gives # levels at each cell?
depth = zeros(size(mbathy));
for z=1:75
    id = (mbathy==z);
    depth(id) = double(nav_lev(z));
end

%% regrid from IPSL grid to CMIP6 grid

figure
pcolor(LAT)
shading flat
figure
pcolor(LON)
shading flat
figure
pcolor(nav_lat)
shading flat
figure
pcolor(nav_lon)
shading flat

%%
dep = griddata(double(nav_lat),double(nav_lon),depth,LAT,LON);

% vdep=depth(:);
% vlon=double(nav_lon(:));
% vlat=double(nav_lat(:));
% dep = griddata(vlon,vlat,vdep,LON,LAT);

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(double(nav_lat),double(nav_lon),depth)
title('IPSL')

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,dep)
title('CMIP')

%%
dep(dep==0) = NaN;
DID = find(~isnan(dep(:))); 
IDN = length(DID); %43942

eq2 = (WID==DID);
sum(eq2)

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = dep(ID);
GRD.lmask = lmask(ID);

GRD.Z(isnan(GRD.Z)) = 0;

%% Save needed variables
save([Cdir 'IPSL/gridspec_ipsl_cm6a.mat'],'depth','nav_lat','nav_lon');
save([Cdir 'IPSL/gridspec_ipsl_cmip6.mat'],'dep','LAT','LON','lmask');
save([Cdir 'IPSL/Data_grid_ipsl.mat'],'GRD');
