% Read Fish-MIP RCP 8.5 netcdfs
% Concatenate bottom temp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CESM/RCP85/';

%%
ncid = netcdf.open([fpath 'cesm_rcp85_o2_zall_monthly_209001-209912.nc4'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 9.96921e+36) = NaN;']);
end
netcdf.close(ncid);

