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

%% Depth
ncid = netcdf.open([Cdir 'cesm_depth.nc4'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%% Retain only water cells
ID = find(MASK(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = H(ID);
GRD.AREA  = AREA(ID);
GRD.lmask = MASK(ID);

%%
ncid = netcdf.open([Cdir 'wc15n_dist.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nj] = size(LON);
DIST = dist(2:(end-1),2:(end-1));
GRD.DIST = DIST(ID);

%% Save needed variables
save([Cdir 'gridspec_3km.mat'],'AREA','H','LAT','LON','MASK','DIST');
save([Cdir 'Data_grid_3km_hist.mat'],'GRD');
          
        