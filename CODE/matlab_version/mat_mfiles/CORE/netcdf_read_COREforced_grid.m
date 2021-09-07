% Read GFDL CORE-forced netcdfs
% grid static

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%%
ncdisp([fpath 'ocean.static.nc'])

% Dimensions:
% xu_ocean = 360
% yu_ocean = 200
% xt_ocean = 360
% yt_ocean = 200

% area_t
% Size:       360x200
% long_name     = 'tracer cell area'
% units         = 'm^2'

% kmt
% Size:       360x200
% long_name     = 'number of depth levels on t-grid'
% _FillValue    = -1.000000020040877e+20
% 
% kmu
% Size:       360x200
% long_name     = 'number of depth levels on u-grid'
% _FillValue    = -1.000000020040877e+20
% 
% ht
% Size:       360x200
% Dimensions: xt_ocean,yt_ocean
% long_name     = 'ocean depth on t-cells'
% units         = 'm'
% _FillValue    = -1.000000020040877e+20
% coordinates   = 'geolon_t geolat_t'
% standard_name = 'sea_floor_depth_below_geoid'
% 
% hu
% Size:       360x200
% Dimensions: xu_ocean,yu_ocean
% long_name     = 'ocean depth on u-cells'
% units         = 'm'
% _FillValue    = -1.000000020040877e+20
% coordinates   = 'geolon_c geolat_c'


%%
ncid = netcdf.open([fpath 'ocean.static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

%%
ht(ht <= -9e19) = NaN;
hu(hu <= -9e19) = NaN;
kmt(kmt <= -9e19) = NaN;
kmu(kmu <= -9e19) = NaN;

%%
figure
pcolor(geolon_t,geolat_t,ht)
shading flat
colorbar

figure
pcolor(geolon_c,geolat_c,hu)
shading flat
colorbar

figure
pcolor(geolon_t,geolat_t,kmt)
shading flat
colorbar

figure
pcolor(geolon_c,geolat_c,kmu)
shading flat
colorbar

%%
save([fpath 'ocean_cobalt_grid.mat'],'geolat_t','geolon_t','ht','kmt','area_t');

