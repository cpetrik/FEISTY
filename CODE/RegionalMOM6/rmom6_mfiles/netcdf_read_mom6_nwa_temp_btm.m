% Read GFDL MOM6-NWA12
% btm temp

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';

%% one file
ncdisp([fpath 'btm_temp.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'])

%%
% Dimensions:
% time = 324   (UNLIMITED)
% yh   = 845
% xh   = 775
% nv   = 2

% btm_temp
% Size:       775x845x324
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
btm_temp_units     = 'deg C';
btm_temp_long_name = 'Bottom Temperature';
FillValue     = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;

% time
time_units    = 'days since 1993-01-01 00:00:00';
calendar_type = 'GREGORIAN';
calendar      = 'gregorian';

%%
ncid = netcdf.open([fpath 'btm_temp.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end

netcdf.close(ncid);

%%
btm_temp(btm_temp>1e19) = nan;
btm_temp = double(btm_temp);

%% Time
yr = 1993 + (time/365);

%%
save([fpath 'btm_temp.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'], '-v7.3');

