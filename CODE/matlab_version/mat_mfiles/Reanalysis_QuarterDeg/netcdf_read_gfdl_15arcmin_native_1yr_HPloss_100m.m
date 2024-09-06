% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss integrated top 100m
% 1/4 grid from Xiao

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';


%% 
ncdisp([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'])

%%
% Size:       1440x1080x600
% Dimensions: xh,yh,z_l,time

% Time
% Size:       12x1
% long_name      = 'time'
time_units       = 'days since 1959-01-01 00:00:00';
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'

% xh
% Size:       1440x1
% Dimensions: xh
% Datatype:   double
% Attributes:
xh_long_name      = 'h point nominal longitude';
xh_units          = 'degrees_east';
% cartesian_axis = 'X'

% yh
% Size:       1080x1
% Dimensions: yh
% Datatype:   double
% Attributes:
yh_long_name      = 'h point nominal latitude';
yh_units          = 'degrees_north';
% cartesian_axis = 'Y'

% jhploss_nmdz_100
% Size:       1440x1080x600
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
nmdz_long_name = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m';
nmdz_units     = 'mol m-2 s-1';
missing_value  = 1.000000020040877e+20;
FillValue      = 1.000000020040877e+20;
% cell_methods  = 'area:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'

% jhploss_nlgz_100
% Size:       1440x1080x600
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
nlgz_long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m';
nlgz_units         = 'mol m-2 s-1';


%% all but HPloss
ncid = netcdf.open([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 3:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% HPloss
for n = 1
    varname = netcdf.inqVar(ncid, n-1);
    jhploss_nlgz_100 = netcdf.getVar(ncid,n-1);
    %jhploss_nlgz_100 = netcdf.getVar(ncid,n-1,[0,0,run1(1)-1],[1440 1080 length(run1)]);
end
for n = 2
    varname = netcdf.inqVar(ncid, n-1);
    jhploss_nmdz_100 = netcdf.getVar(ncid,n-1);
    %jhploss_nmdz_100 = netcdf.getVar(ncid,n-1,[0,0,run1(1)-1],[1440 1080 length(run1)]);
end
netcdf.close(ncid);

jhploss_nmdz_100(jhploss_nmdz_100 >= 1.00e+19) = NaN;
jhploss_nlgz_100(jhploss_nlgz_100 >= 1.00e+19) = NaN;

inan = isnan(jhploss_nmdz_100(:));
ocell = ~isnan(squeeze(jhploss_nmdz_100(:,:,1,1)));
lcell = isnan(squeeze(jhploss_nmdz_100(:,:,1,1)));

%% Time
yr = 1959 + (time/365);

%% grid - original
ncid = netcdf.open([spath 'full.ocean_static_FishMIP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
test = double(squeeze(jhploss_nmdz_100(:,:,6)));
test2 = double(squeeze(jhploss_nlgz_100(:,:,6)));

%% map
clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%%
figure
pcolor(geolon,geolat,test); shading flat; colorbar; 
colormap('jet')

figure
pcolor(geolon,geolat,test2); shading flat; colorbar;
colormap('jet')

geolon = double(geolon);
geolat = double(geolat);
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',mlonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,test)
colormap('jet')
colorbar
%clim([0 0.01])

figure
pcolor(test2); shading flat; colorbar;
colormap('jet')

% figure
% pcolor(LON,LAT,test); shading flat; colorbar;
% colormap('jet')
% 
% figure
% pcolor(LON,LAT,test2); shading flat; colorbar;
% colormap('jet')


