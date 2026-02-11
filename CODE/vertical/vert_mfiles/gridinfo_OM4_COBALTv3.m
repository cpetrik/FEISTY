% 0.5 degree global MOM6-COBALTv3 only


clear
close all

fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.static.nc'])

areacello_units      = 'm2';
areacello_long_name  = 'Ocean Grid-Cell Area';
deptho_units      = 'm';
deptho_long_name  = 'Sea Floor Depth';
dxt_units         = 'm';
dxt_long_name     = 'Delta(x) at thickness/tracer points (meter)';
dyt_units         = 'm';
dyt_long_name     = 'Delta(y) at thickness/tracer points (meter)';
geolat_units      = 'degrees_north';
geolat_long_name  = 'Latitude of tracer (T) points';
geolon_units      = 'degrees_east';
geolon_long_name  = 'Longitude of tracer (T) points';
sftof_units       = '%';
sftof_long_name   = 'Sea Area Fraction';
wet_long_name     = '0 if land, 1 if ocean at tracer points';

%%
ncdisp([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.volcello.nc'])

z_l_units     = 'meters';
z_l_long_name = 'Depth at cell center';

%% static lat lon
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

%% vertical
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.volcello.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for n = 1 %:nvars
    varname = netcdf.inqVar(ncid, n-1);
    z_l = netcdf.getVar(ncid,n-1);
end
netcdf.close(ncid);

%% remove fill value
deptho(deptho>1e19) = nan;

%% 
save([fpath 'grid_OM4_05_COBALTv3.mat'],'areacello_units','areacello_long_name','areacello',...
    'deptho_units','deptho_long_name','deptho','dxt_units','dxt_long_name','dxt',...
    'dyt_units','dyt_long_name','dyt','geolat_units','geolat_long_name','geolat',...
    'geolon_units','geolon_long_name','geolon','sftof_units','sftof_long_name','sftof',...
    'wet_long_name','wet','z_l_units','z_l_long_name','z_l');










