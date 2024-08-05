% Read GFDL 1/4 netcdfs
% nitrogen prim. prod. integrals in upper 100m
% Raw units provided by GFDL
% Before conversion by Matthias

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_cobalt_fluxes_int_FishMIP.nc'])

%%
% Global Attributes:
% filename         = '20000101.ocean_cobalt_fluxes_int.nc'
% title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% associated_files = 'areacello: 20000101.ocean_static.nc'
% grid_type        = 'regular'
% grid_tile        = 'N/A'
% history          = 'Thu Mar 10 09:35:08 2022: ncks -v jprod_nlgp_100,jprod_ndi_100,jprod_nsmp_100 20000101.ocean_cobalt_fluxes_int.nc 20000101.ocean_cobalt_fluxes_int_FishMIP.nc'
% NCO              = '"4.5.4"'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20

% Dimensions:
% time = 12    (UNLIMITED)
% yh   = 1080
% xh   = 1440

% jprod_ndi_100
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'Diazotroph nitrogen prim. prod. integral in upper 100m'
% units         = 'mol m-2 s-1'

% jprod_nlgp_100
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'Large phyto. nitrogen  prim. prod. integral in upper 100m'
% units         = 'mol m-2 s-1'

% jprod_nsmp_100
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'Small phyto. nitrogen  prim. prod. integral in upper 100m'
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

%%
ncid = netcdf.open([fpath '20000101.ocean_cobalt_fluxes_int_FishMIP.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);


%% Time
yr = ((time+1)/12)+1601;


%%
test=(double(squeeze(jprod_ndi_100(:,:,6))));

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
