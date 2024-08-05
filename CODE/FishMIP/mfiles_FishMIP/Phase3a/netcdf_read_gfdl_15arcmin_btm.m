% Read GFDL 1/4 netcdfs
% Bottom temp and ndet sinking flux to bottom
% Raw units provided by GFDL
% Before conversion by Matthias

clear 
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_cobalt_btm_FishMIP.nc'])

%%
% Global Attributes:
% filename         = '20000101.ocean_cobalt_btm.nc'
% title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% associated_files = 'areacello: 20000101.ocean_static.nc'
% grid_type        = 'regular'
% grid_tile        = 'N/A'
% history          = 'Thu Mar 10 09:35:12 2022: ncks -v fndet_btm,btm_temp 20000101.ocean_cobalt_btm.nc 20000101.ocean_cobalt_btm_FishMIP.nc'
% NCO              = '"4.5.4"'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20

% Dimensions:
% time = 12    (UNLIMITED)
% yh   = 1080
% xh   = 1440

% btm_temp
% Size:       1440x1080x12
% long_name     = 'Bottom Temperature'
% units         = 'deg C'

% fndet_btm
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'ndet sinking flux to bottom'
% units         = 'mol m-2 s-1'

% time
% Size:       12x1
% Dimensions: time
% long_name      = 'time'
% units          = 'days since 1959-01-01 00:00:00'
% cartesian_axis = 'T'
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'
% bounds         = 'time_bnds'

% xh
% Size:       1440x1
% Dimensions: xh
% long_name      = 'h point nominal longitude'
% units          = 'degrees_east'
% cartesian_axis = 'X'

% yh
% Size:       1080x1
% Dimensions: yh
% Datatype:   double
% Attributes:
% long_name      = 'h point nominal latitude'
% units          = 'degrees_north'
% cartesian_axis = 'Y'

% z_l
% Size:       35x1
% Dimensions: z_l
% Datatype:   double
% Attributes:
% long_name      = 'Depth at cell center'
% units          = 'meters'

%%
ncid = netcdf.open([fpath '20000101.ocean_cobalt_btm_FishMIP.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);


%%
test=(double(squeeze(fndet_btm(:,:,6))));

figure
pcolor(test); shading flat;
title('pcolor(test)')

%%
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
ncid = netcdf.open([cpath 'full.ocean_static_FishMIP.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

geolat = double(geolat);
geolon = double(geolon);
figure
pcolor(geolat,geolon,test); shading flat;
title('pcolor(geolat,geolon,test)')

