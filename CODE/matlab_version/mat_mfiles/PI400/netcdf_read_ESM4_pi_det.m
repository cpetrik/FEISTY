% Read GFDL ESM4 Preindust netcdfs
% bottom POC flux

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/ESM4_PI/';

%%
ncdisp([fpath 'ocean_cobalt_btm.100101-102012.fndet_btm.nc'])

% Dimensions:
% xu_ocean = 360
% yu_ocean = 200
% time     = 240   (UNLIMITED)
% nv       = 2
% xt_ocean = 360
% yt_ocean = 200
% time:units = 'days since 0001-01-01 00:00:00'

% fndet_btm
% Size:       360x200x240
% Dimensions: xt_ocean,yt_ocean,time
% Datatype:   single
% Attributes:
long_name     = 'ndet sinking flux to bottom';
units         = 'mol m-2 s-1';
% missing_value = -10000000000
% _FillValue    = -10000000000
% cell_methods  = 'time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% coordinates   = 'geolon_t geolat_t'

%% 1st 20 yrs
ncid = netcdf.open([fpath 'ocean_cobalt_btm.100101-102012.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

%% NaNs on land cells
fndet_btm(fndet_btm <= -1e9) = NaN;
