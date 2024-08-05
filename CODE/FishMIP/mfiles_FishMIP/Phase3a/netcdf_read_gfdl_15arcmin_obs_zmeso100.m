% Read GFDL netcdfs
% obsclim
% zmeso all depths, integrate top 100m
% 1/4 grid from ISIMIP

clear 
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

%% zmeso
ncdisp([fpath 'gfdl-mom6-cobalt2_obsclim_zmeso_15arcmin_global_monthly_1961_2010.nc'])

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
% 
% lev
% Size:       35x1
% Dimensions: lev
% Datatype:   double
% Attributes:
lev_long_name    = 'Depth at cell center';
lev_units        = 'meters';
% positive       = 'down'
% axis           = 'Z'
% cartesian_axis = 'Z'
% edges          = 'z_i'
% 
% zmeso
% Size:       1440x720x35x600
% Dimensions: lon,lat,lev,time
% Datatype:   single
% Attributes:
standard_name        = 'mole_concentration_of_mesoplankton_expressed_as_carbon_in_sea_water';
long_name            = 'Mesozooplankton Carbon Concentration';
units                = 'mol m-3';
FillValue            = 1.000000020040877e+20;
missing_value        = 1.000000020040877e+20;
GFDL_variable        = 'nlgz nmdz';
aggregated_variables = 'zoolarge zoomedium';

%% Cell thickness = thkcello
tcid = netcdf.open([fpath 'gfdl-mom6-cobalt2_obsclim_thkcello_15arcmin_global_fixed.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    eval([ varname ' = netcdf.getVar(tcid,t-1);']);
    %thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;
thkcello = double(thkcello);

%% all but zmeso
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_obsclim_zmeso_15arcmin_global_monthly_1961_2010.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subset of depth
z100 = find(lev <= 100);
thkcello = thkcello(:,:,z100);

% do one time half and then the other
run1 = 1:300;
run2 = 301:600;

%% zmeso 1st half
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1,[0,0,0,run1(1)-1],[1440 720 length(z100) length(run1)]);
end
zmeso(zmeso >= 1.00e+19) = NaN;

%% Integrate top 100 m
thk3 = repmat(nansum(thkcello,3),1,1,length(run1));
zmeso_100 = squeeze(nansum((zmeso.*thkcello),3));

%% zmeso 2nd half
clear zmeso

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1,[0,0,0,run2(1)-1],[1440 720 length(z100) length(run2)]);
end
netcdf.close(ncid);

%% Integrate top 100 m
zmeso2_100 = squeeze(nansum((zmeso.*thkcello),3)) ;

%% grid
[LAT,LON] = meshgrid(lat, lon);

%% Time
yr = 1901 + (time/12);

%%
clear zmeso
zmeso_100(:,:,run2) = zmeso2_100;
zmeso_100 = double(zmeso_100);

mtp = squeeze(nanmean(nanmean(zmeso_100,2),1));
plot(yr,mtp)

%%
clear zmeso2_100
%%
units_orig = units;
units_vint = 'mol m-2';

%%
save([fpath 'gfdl-mom6-cobalt2_obsclim_zmeso100_15arcmin_global_monthly_1961_2010.mat'],...
    'GFDL_variable','long_name','standard_name','missing_value','units',...
    'lat','lon','time','LAT','LON','zmeso_100','time_units','yr',...
    'run1','run2','z100','lev','lev_long_name','lev_units','-v7.3');



