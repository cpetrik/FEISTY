% Read GFDL 1/4 netcdfs
% obsclim
% Bottom temp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/';
%fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/';

%% one file
ncdisp([fpath 'gfdl-mom6-cobalt2_obsclim_tob_15arcmin_global_monthly_1961_2010.nc'])

%%
% time
% Size:       600x1
% Dimensions: time
% Datatype:   double
% Attributes:
% standard_name = 'time'
% long_name     = 'time'
time_units      = 'months since 1901-1-1 00:00:00';
% calendar      = '360_day'
% axis          = 'T'

% tob
% Size:       1440x720x600
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature on Bottom';
units         = 'degC';
FillValue     = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
GFDL_variable = 'btm_temp';

%%
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_obsclim_tob_15arcmin_global_monthly_1961_2010.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(tcid);

tob(tob>1e19) = nan;
tob = double(tob);

%% grid
[LAT,LON] = meshgrid(lat, lon);

%% Time
yr = 1901 + (time/12);

%%
save([fpath 'QuarterDeg/gfdl-mom6-cobalt2_obsclim_tob_15arcmin_global_monthly_1961_2010.mat'],...
    'GFDL_variable','long_name','standard_name','missing_value','units',...
    'lat','lon','time','LAT','LON','tob','time_units','yr','-v7.3');
