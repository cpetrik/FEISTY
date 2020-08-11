% Make GRD file for FEISTY input from CESM 1 degree model

clear all
close all

Cdir = '/Volumes/FEISTY/Fish-MIP/CMIP6/';

%% Lat & Lon
ncid = netcdf.open([Cdir 'IPSL/hist/ipsl-cm6a-lr_r1i1p1f1_historical_tob_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
tob(tob >= 1.00e+20) = NaN;

%Land mask
mask = squeeze(tob(:,:,1,1));
lmask = mask;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);
%Lon=0 is prime meridian

%% Depth
ncid = netcdf.open([Cdir 'GFDL/ocean_cobalt_omip_tracers_month_z_1x1deg.static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
deptho(deptho >= 1.00e+20) = NaN;

%% Retain only water cells
ID = find(lmask(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = double(LON(ID));
GRD.LAT = double(LAT(ID));
GRD.Z   = double(depth(ID));
GRD.lmask = double(lmask(ID));
GRD.AREA  = double(areacello(ID));

%% Save needed variables
save([Cdir 'IPSL/gridspec_ipsl.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'IPSL/Data_grid_ipsl.mat'],'GRD');
