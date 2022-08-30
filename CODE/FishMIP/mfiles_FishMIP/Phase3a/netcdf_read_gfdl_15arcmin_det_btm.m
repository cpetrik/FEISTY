% Read GFDL 1/4 netcdfs
% ctrlclim
% Detritus sinking flux btm

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/';
fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/';

%% one file
ncdisp([fpath 'gfdl-mom6-cobalt2_ctrlclim_expc-bot_15arcmin_global_monthly_1961_2010.nc'])

%%
% time
% Size:       600x1
% Dimensions: time
% Datatype:   double
% Attributes:
% standard_name = 'time'
% long_name     = 'time'
% units         = 'months since 1901-1-1 00:00:00'
% calendar      = '360_day'
% axis          = 'T'
% 
% expc-bot
% Size:       1440x720x600
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
long_name     = 'Sinking Particulate Organic Carbon Flux on Bottom (z_b)';
units         = 'mol m-2 s-1';
FillValue    = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
GFDL_variable = 'fndet_btm';

%%
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_ctrlclim_expc-bot_15arcmin_global_monthly_1961_2010.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    det_btm = netcdf.getVar(ncid,n-1);
    
end
netcdf.close(ncid);

%%
det_btm(det_btm>1e19) = nan;
det_btm = double(det_btm);

%% grid
[LAT,LON] = meshgrid(lat, lon);

%% Time
yr = 1901 + (time/12);

%%
save([fpath 'QuarterDeg/gfdl-mom6-cobalt2_ctrlclim_expc-bot_15arcmin_global_monthly_1961_2010.mat'],...
    'GFDL_variable','long_name','standard_name','missing_value','units',...
    'lat','lon','time','LAT','LON','det_btm','yr', '-v7.3');

