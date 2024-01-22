% Make GRD file for FEISTY input from CESM 1 degree model

clear all
close all

Cdir = '/Volumes/FEISTY/Fish-MIP/CMIP6/';

%% Lat & Lon
ncid = netcdf.open([Cdir 'GFDL/hist/gfdl-esm4_r1i1p1f1_historical_thetao-bot_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
thetao(thetao >= 1.00e+20) = NaN;

%Land mask
mask = squeeze(thetao(:,:,1));
lmask = double(mask);
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);

LID = find(lmask(:)==1);
WID = find(~isnan(thetao(:,:,1)));  % spatial index of water cells
NID = length(WID);

eq1 = (WID==LID);
sum(eq1)

%% Depth, lat, lon
ncid = netcdf.open([Cdir 'GFDL/ocean_cobalt_omip_tracers_month_z_1x1deg.static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
deptho(deptho >= 1.00e+20) = NaN;
geolon(geolon >= 1.00e+20) = NaN;
geolat(geolat >= 1.00e+20) = NaN;
wet(wet >= 1.00e+20) = 0;

%% regrid (shift) from CM4 grid to CMIP6 grid
deptho = double(deptho);
areac = double(areacello);

geolon_orig = geolon; 
test = geolon; 
id=find(test < -180);
test(id)=test(id)+360;
geolon = double(test);

%%
mlon=find(geolon(:,end)<-179);

depth(1:180,:) = deptho(mlon:end,:);
depth(181:360,:) = deptho(1:(mlon-1),:);

area(1:180,:) = areac(mlon:end,:);
area(181:360,:) = areac(1:(mlon-1),:);

DID = find(~isnan(depth(:))); 
IDN = length(DID);

eq2 = (WID==DID);
sum(eq2)

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = depth(ID);
GRD.lmask = lmask(ID);
GRD.AREA  = area(ID);

%% Save needed variables
save([Cdir 'GFDL/gridspec_gfdl_cm4.mat'],'deptho','geolat','geolon','wet','areac');
save([Cdir 'GFDL/gridspec_gfdl_cmip6.mat'],'depth','LAT','LON','lmask','area');
save([Cdir 'GFDL/Data_grid_gfdl.mat'],'GRD');
