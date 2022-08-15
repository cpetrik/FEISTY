% Read GFDL 1/4 netcdfs
% ctrlclim
% temp all depths
% take mean over top 100m

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
%fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/';

%% one file
ncdisp([fpath 'gfdl-mom6-cobalt2_ctrlclim_thetao_onedeg_global_monthly_1961_2010.nc'])

%%
% time
% Size:       600x1
% Dimensions: time
% Datatype:   double
% Attributes:
% standard_name = 'time'
% long_name     = 'time'
time_units         = 'months since 1901-1-1 00:00:00';
% calendar      = '360_day'
% axis          = 'T'

% Dimensions: lev
% Datatype:   double
% Attributes:
lev_long_name      = 'Depth at cell center';
lev_units          = 'meters';
% positive       = 'down'
% axis           = 'Z'
% cartesian_axis = 'Z'
% edges          = 'z_i'

% thetao
% Size:       1440x720x35x600
% Dimensions: lon,lat,lev,time
% Datatype:   single
% Attributes:
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';
FillValue     = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
GFDL_variable = 'thetao';

%% Cell thickness = thkcello
tcid = netcdf.open([fpath 'gfdl-mom6-cobalt2_ctrlclim_thkcello_onedeg_global_fixed.nc'],'NC_NOWRITE');
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

%% all but temp
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_ctrlclim_thetao_onedeg_global_monthly_1961_2010.nc'],'NC_NOWRITE');
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

%% temp 1st half
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    thetao = netcdf.getVar(ncid,n-1,[0,0,0,run1(1)-1],[360 180 length(z100) length(run1)]);
end
thetao(thetao >= 1.00e+19) = NaN;

%% Mean top 100 m
thk3 = repmat(nansum(thkcello,3),1,1,length(run1));
temp_100 = squeeze(nansum((thetao.*thkcello),3)) ./ thk3;

%% temp 2nd half
clear thetao

ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_ctrlclim_thetao_onedeg_global_monthly_1961_2010.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    thetao = netcdf.getVar(ncid,n-1,[0,0,0,run2(1)-1],[360 180 length(z100) length(run2)]);
end
netcdf.close(ncid);

%% Mean top 100 m
temp2_100 = squeeze(nansum((thetao.*thkcello),3)) ./ thk3;

%% grid
[LAT,LON] = meshgrid(lat, lon);

%% Time
yr = 1901 + (time/12);

%%
clear thetao
temp_100(:,:,run2) = temp2_100;
temp_100 = double(temp_100);

mtp = squeeze(nanmean(nanmean(temp_100,2),1));
plot(yr,mtp)

%%
clear temp2_100

save([fpath 'gfdl-mom6-cobalt2_ctrlclim_temp100_onedeg_global_monthly_1961_2010.mat'],...
    'GFDL_variable','long_name','standard_name','missing_value','units',...
    'lat','lon','time','LAT','LON','temp_100','time_units','yr',...
    'run1','run2','z100','lev','lev_long_name','lev_units','-v7.3');
