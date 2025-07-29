% Read GFDL MOM6-NWA12
% zmeso integrated over top 100m

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';

%% zmeso
ncdisp([fpath 'jhploss_nmdz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'])

ncdisp([fpath 'jhploss_nlgz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'])

%%
% Dimensions:
% time = 324   (UNLIMITED)
% yh   = 845
% xh   = 775
% nv   = 2

% jhploss_nmdz_100
% Size:       775x845x324
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
FillValue    = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
jhploss_nmdz_units     = 'mol m-2 s-1';
jhploss_nmdz_long_name = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m';


% jhploss_nlgz_100
jhploss_nlgz_units     = 'mol m-2 s-1';
jhploss_nlgz_long_name = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m';

% time
time_units    = 'days since 1993-01-01 00:00:00';
calendar_type = 'GREGORIAN';
calendar      = 'gregorian';

%% medium
ncid = netcdf.open([fpath 'jhploss_nmdz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% large
ncid = netcdf.open([fpath 'jhploss_nlgz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

netcdf.close(ncid);

%%
jhploss_nmdz_100(jhploss_nmdz_100>1e19) = nan;
jhploss_nlgz_100(jhploss_nlgz_100>1e19) = nan;

jhploss_nmdz_100 = double(jhploss_nmdz_100);
jhploss_nlgz_100 = double(jhploss_nlgz_100);

%% Time
yr = 1993 + (time/365);

%%
save([fpath 'jhploss_mdz_lgz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],...
    '-v7.3');



