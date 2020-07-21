% Read Fish-MIP netcdfs
% 

clear all
close all

%%
fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
ncid = netcdf.open([fpath 'cesm_pi_lzoo_zall_monthly_185001-185912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
ncid = netcdf.open([fpath 'cesm_pi_to_zall_monthly_185001-185912.nc4'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
ncid = netcdf.open([fpath 'cesm_pi_to_zb_monthly_185001-185912.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%%
fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
ncid = netcdf.open([fpath 'cesm_pi_POC_FLUX_IN_zall_monthly_185001-185912.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NO BOTTOM POC FLUX, WOULD NEED TO GO THROUGH EACH LAYER FROM BOTTOM TO TOP
%TO FILL IN

%% "POC Production - POC_PROD:units = "mmol/m^3/s"
fpath='/Volumes/GFDL/Fish-MIP/CESM/PreIndust/';
ncid = netcdf.open([fpath 'cesm_pi_POC_PROD_zall_monthly_185001-185912.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

%NO BOTTOM POC FLUX, WOULD NEED TO GO THROUGH EACH LAYER FROM BOTTOM TO TOP
%TO FILL IN

%% Calculate quantities
dat = dxt.*dyt;
datr = 1.0./(dat+eps);

tmask = wet;

%% Save needed variables
save([fpath fname '.mat'])

%% Retain all 2D
GRD.LON = double(geolon);
GRD.LAT = double(geolat);
GRD.Z   = double(deptho);
GRD.DX = double(dxt);
GRD.DY = double(dyt);
GRD.AREA  = double(areacello);
GRD.dxtn  = double(dxt);
GRD.dyte  = double(dyt);
GRD.datr  = double(datr);
GRD.lmask = double(tmask);

%! save
save([fpath 'Data_grid2D_MOM6_preindust.mat'],'GRD');

%% Retain 2D for mom6 diffusion
[ni,nj] = size(dat);
G.isc = 1;
G.iec = ni;
G.jsc = 2; %1st grid cell is Antarctica
G.jec = nj;
%don't know what these are, assume same as above
G.IscB = 1;
G.IecB = ni;
G.JscB = 2;
G.JecB = nj;
G.dxCu = double(dxCu);
G.dxCv = double(dxCv);
G.dyCu = double(dyCu);
G.dyCv = double(dyCv);
G.area = double(areacello);
G.Iarea = double(datr);
G.dxt  = double(dxt);
G.dyt  = double(dyt);
G.mask = double(tmask);

%! save
save([fpath 'Data_grid2D_diff_MOM6_preindust.mat'],'G');

%% Retain only water cells and vectorize
clear GRD
ID = find(wet(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = geolon(ID);
GRD.LAT = geolat(ID);
GRD.Z   = deptho(ID);
GRD.DX = dxt(ID);
GRD.DY = dyt(ID);
GRD.AREA  = areacello(ID);
GRD.dxtn  = dxt(ID);
GRD.dyte  = dyt(ID);
GRD.datr  = datr(ID);
GRD.lmask = tmask(ID);

%! save
save([fpath 'Data_grid1D_MOM6_preindust.mat'],'GRD');
        



