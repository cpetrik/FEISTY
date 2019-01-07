% Make GRD file for FEISTY input

clear all
close all

Cdir = '/Volumes/GFDL/NEMURO/3km/';

ncid = netcdf.open([Cdir 'wc15n_avgmo_1988.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 99999) = NaN;']);
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
          
        