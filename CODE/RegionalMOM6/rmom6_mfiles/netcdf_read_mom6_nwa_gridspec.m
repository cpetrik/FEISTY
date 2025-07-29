% Read GFDL 1/4 netcdfs
% ctrl thkcello

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';

%%
ncdisp([fpath 'ocean_static.nc'])

%%
FillValue    = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;

% areacello
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
areacello_units         = 'm2';
areacello_long_name     = 'Ocean Grid-Cell Area';


% areacello_cu
% Size:       776x845
% Dimensions: xq,yh
% Datatype:   single
areacello_cu_units         = 'm2';
areacello_cu_long_name     = 'Ocean Grid-Cell Area at u points';

% areacello_cv
% Size:       775x846
% Dimensions: xh,yq
% Datatype:   single
areacello_cv_units         = 'm2';
areacello_cv_long_name     = 'Ocean Grid-Cell Area at v points';

% deptho
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
deptho_units         = 'm';
deptho_long_name     = 'Sea Floor Depth';
deptho_standard_name = 'sea_floor_depth_below_geoid';

% dxCu
% Size:       776x845
% Dimensions: xq,yh
% Datatype:   single
dxCu_units         = 'm';
dxCu_long_name     = 'Delta(x) at u points (meter)';

% dxCv
% Size:       775x846
% Dimensions: xh,yq
% Datatype:   single
units         = 'm';
long_name     = 'Delta(x) at v points (meter)';

% dxt
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
dxCv_units         = 'm';
dxCv_long_name     = 'Delta(x) at thickness/tracer points (meter)';

% dyCu
% Size:       776x845
% Dimensions: xq,yh
% Datatype:   single
dyCu_units         = 'm';
dyCu_long_name     = 'Delta(y) at u points (meter)';

% dyCv
% Size:       775x846
% Dimensions: xh,yq
% Datatype:   single
dyCv_units         = 'm';
dyCv_long_name     = 'Delta(y) at v points (meter)';

% dyt
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
dyt_units         = 'm';
dyt_long_name     = 'Delta(y) at thickness/tracer points (meter)';

% geolat
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
geolat_units         = 'degrees_north';
geolat_long_name     = 'Latitude of tracer (T) points';

% geolat_u
% Size:       776x845
% Dimensions: xq,yh
% Datatype:   single
geolat_u_units         = 'degrees_north';
geolat_u_long_name     = 'Latitude of zonal velocity (Cu) points';

% geolat_v
% Size:       775x846
% Dimensions: xh,yq
% Datatype:   single
geolat_v_units         = 'degrees_north';
geolat_v_long_name     = 'Latitude of meridional velocity (Cv) points';

% geolon
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
geolon_units         = 'degrees_east';
geolon_long_name     = 'Longitude of tracer (T) points';

% geolon_u
% Size:       776x845
% Dimensions: xq,yh
% Datatype:   single
geolon_u_units         = 'degrees_east';
geolon_u_long_name     = 'Longitude of zonal velocity (Cu) points';

% geolon_v
% Size:       775x846
% Dimensions: xh,yq
% Datatype:   single
geolon_v_units         = 'degrees_east';
geolon_v_long_name     = 'Longitude of meridional velocity (Cv) points';

% sftof
% Size:       775x845
% Dimensions: xh,yh
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
% Size:       775x845
% Dimensions: xh,yh
% Datatype:   single
wet_long_name     = '0 if land, 1 if ocean at tracer points';

% wet_u
% Size:       776x845
% Dimensions: xq,yh
% Datatype:   single
wet_u_long_name     = '0 if land, 1 if ocean at zonal velocity (Cu) points';

% wet_v
% Size:       775x846
% Dimensions: xh,yq
% Datatype:   single
wet_v_long_name     = '0 if land, 1 if ocean at meridional velocity (Cv) points';

% xh
% Size:       775x1
% Dimensions: xh
% Datatype:   double
xh_units     = 'degrees_east';
xh_long_name = 'h point nominal longitude';
xh_axis      = 'X';

% xq
% Size:       776x1
% Dimensions: xq
% Datatype:   double
xq_units     = 'degrees_east';
xq_long_name = 'q point nominal longitude';
xq_axis      = 'X';

% yh
% Size:       845x1
% Dimensions: yh
% Datatype:   double
yh_units     = 'degrees_north';
yh_long_name = 'h point nominal latitude';
yh_axis      = 'Y';

% yq
% Size:       846x1
% Dimensions: yq
% Datatype:   double
yq_units     = 'degrees_north';
yq_long_name = 'q point nominal latitude';
yq_axis      = 'Y';

%%
ncid = netcdf.open([fpath 'ocean_static.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%%
save([fpath 'nwa_raw_ocean_static_gridspec.mat']);



