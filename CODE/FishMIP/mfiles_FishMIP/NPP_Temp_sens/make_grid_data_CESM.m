% Make GRD file for FEISTY input from CESM 1 degree model

clear all
close all

Cdir = '/Volumes/GFDL/Fish-MIP/CESM/';

%% Lat & Lon
ncid = netcdf.open([Cdir 'Hist/cesm_hist_szoo_zall_monthly_185001-185912.nc4'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
szoo(szoo>=9e+36)=nan;

%Land mask
mask = squeeze(szoo(:,:,1,1));
lmask = mask;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);
%Lon=0 is prime meridian

%% Depth
ncid = netcdf.open([Cdir 'cesm_depth.nc4'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%depth in cm, change to m
depth = HT*1e-2;
depth(depth>=9e+36)=nan;

%% Retain only water cells
ID = find(lmask(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = depth(ID);
GRD.lmask = lmask(ID);

%% Save needed variables
save([Cdir 'gridspec_cesm.mat'],'depth','LAT','LON','lmask','mask');
save([Cdir 'Data_grid_cesm.mat'],'GRD');
          
        