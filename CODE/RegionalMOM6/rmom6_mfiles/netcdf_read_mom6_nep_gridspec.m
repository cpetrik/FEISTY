% Read GFDL 1/4 netcdfs
% ctrl thkcello

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';

%%
ncdisp([fpath 'ocean_static.nc'])

%%
FillValue    = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;

% areacello
% Size:       343x817
% Dimensions: ih,ih
% Datatype:   single
areacello_units         = 'm2';
areacello_long_name     = 'Ocean Grid-Cell Area';

% deptho
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
deptho_units         = 'm';
deptho_long_name     = 'Sea Floor Depth';
deptho_standard_name = 'sea_floor_depth_below_geoid';

% dxt
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
dxCv_units         = 'm';
dxCv_long_name     = 'Delta(x) at thickness/tracer points (meter)';

% dyt
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
dyt_units         = 'm';
dyt_long_name     = 'Delta(y) at thickness/tracer points (meter)';

% geolat
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
geolat_units         = 'degrees_north';
geolat_long_name     = 'Latitude of tracer (T) points';

% geolat_u
% Size:       342x816
% Dimensions: iq,jh
% Datatype:   single
geolat_u_units         = 'degrees_north';
geolat_u_long_name     = 'Latitude of zonal velocity (Cu) points';

% geolat_v
% Size:       342x816
% Dimensions: ih,jq
% Datatype:   single
geolat_v_units         = 'degrees_north';
geolat_v_long_name     = 'Latitude of meridional velocity (Cv) points';

% geolon
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
geolon_units         = 'degrees_east';
geolon_long_name     = 'Longitude of tracer (T) points';

% geolon_u
% Size:       342x816
% Dimensions: iq,jh
% Datatype:   single
geolon_u_units         = 'degrees_east';
geolon_u_long_name     = 'Longitude of zonal velocity (Cu) points';

% geolon_v
% Size:       342x816
% Dimensions: ih,jq
% Datatype:   single
geolon_v_units         = 'degrees_east';
geolon_v_long_name     = 'Longitude of meridional velocity (Cv) points';

% sftof
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
sftof_units         = '%';
sftof_long_name     = 'Sea Area Fraction';

% time
% Size:       1x1
% Dimensions: time
% Datatype:   double
time_units         = 'days since 1980-01-01 00:00:00';
time_calendar_type = 'GREGORIAN';
time_calendar      = 'gregorian';

% wet
% Size:       342x816
% Dimensions: ih,jh
% Datatype:   single
wet_long_name     = '0 if land, 1 if ocean at tracer points';

% ih
% Size:       342x1
% Dimensions: ih
% Datatype:   double
xh_units     = 'degrees_east';
xh_long_name = 'h point nominal longitude';
xh_axis      = 'X';

% jh
% Size:       816x1
% Dimensions: jh
% Datatype:   double
yh_units     = 'degrees_north';
yh_long_name = 'h point nominal latitude';
yh_axis      = 'Y';


%%
ncid = netcdf.open([fpath 'ocean_static.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%%
save([fpath 'nep_raw_ocean_static_gridspec.mat']);



