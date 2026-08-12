% Read GFDL MOM6-NEP10
% Detritus sinking flux btm

clear
close all

%fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';
fpath='/project/Feisty/GCM_Data/MOM6-NEP10/';

%% one file
ncdisp([fpath 'fntot_btm.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc'])

%%
% Dimensions:
% time = 324   (UNLIMITED)
% yh   = 845
% xh   = 775
% nv   = 2

% fntot_btm
% Size:       775x845x324
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
FillValue     = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
fntot_btm_units      = 'mol m-2 s-1';
fntot_btm_long_name  = 'Total N sinking flux to bottom';
% cell_methods  = 'area:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'

% time
time_units    = 'days since 1993-01-01 00:00:00';
calendar_type = 'GREGORIAN';
calendar      = 'gregorian';

%%
ncid = netcdf.open([fpath 'fntot_btm.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end

netcdf.close(ncid);

%%
fntot_btm(fntot_btm>1e19) = nan;
fntot_btm = double(fntot_btm);

%% Time
yr = 1993 + (time/365);

%%
save([fpath 'fntot_btm.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat'], '-v7.3');

