% Read GFDL netcdfs
% obsclim
% mdz & lgz biom all depths
% 1/4 grid from ISIMIP

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% zoo regular grid
ncdisp([fpath '19630101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc'])

% Global Attributes:
%            CDI              = 'Climate Data Interface version 2.4.2 (https://mpimet.mpg.de/cdi)'
%            Conventions      = 'CF-1.6'
%            filename         = '19620101.ocean_cobalt_tracers_month_z.nc'
%            title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
%            associated_files = 'areacello: 19620101.ocean_static.nc'
%            grid_type        = 'regular'
%            grid_tile        = 'N/A'
%            history          = 'Tue Jul 30 13:08:17 2024: cdo -selname,nlgz,nmdz -sellonlatbox,-180,180,-90,90 -remapbil,grid.des /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19620101.ocean_cobalt_tracers_month_z_FishMIP_CP.nc.tmp /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19620101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc
%                               Tue Jul 30 12:52:30 2024: cdo -O -merge /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19620101.ocean_cobalt_tracers_month_z_FishMIP_CP.nc full.ocean_static_FishMIP.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19620101.ocean_cobalt_tracers_month_z_FishMIP_CP.nc.tmp
%                               Fri Nov 17 12:19:13 2023: ncks -v nmdz,nlgz 19620101.ocean_cobalt_tracers_month_z.nc 19620101.ocean_cobalt_tracers_month_z_FishMIP_CP.nc'
%            NCO              = 'netCDF Operators version 5.0.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)'
%            CDO              = 'Climate Data Operators version 2.4.2 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            time = 12    (UNLIMITED)
%            lon  = 1440
%            lat  = 720
%            z_l  = 35
% Variables:
%     time
%            Size:       12x1
        time_units         = 'days since 1959-01-01 00:00:00';
%
%     lon
%            Size:       1440x1
%            units         = 'degrees_east'
%
%     lat
%            Size:       720x1
%            units         = 'degrees_north'
%
%     z_l
%            Size:       35x1
%            long_name      = 'Depth at cell center'
%                        units          = 'meters'
%
%     nlgz
%            Size:       1440x720x35x12
%            Dimensions: lon,lat,z_l,time
%            Attributes:
        nlgz_long_name     = 'large Zooplankton Nitrogen';
        nlgz_units         = 'mol/kg';
        missing_value = 1.000000020040877e+20;
%
%     nmdz
%            Size:       1440x720x35x12
%            Dimensions: lon,lat,z_l,time
%            Attributes:
        nmdz_long_name     = 'Medium-sized zooplankton Nitrogen';
        nmdz_units         = 'mol/kg';
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'

%%
ncid = netcdf.open([fpath '19630101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = 1959 + (time/365);

%%
nmdz(nmdz >= 1.00e+19) = NaN;
nlgz(nlgz >= 1.00e+19) = NaN;

inan = isnan(nmdz(:));
ocell = ~isnan(squeeze(nmdz(:,:,1,1)));
lcell = isnan(squeeze(nmdz(:,:,1,1)));

%% molN/m3 --> gC/m3
nmdz = double(nmdz * (106/16) * 12.01);
nlgz = double(nlgz * (106/16) * 12.01);

%% grid - remapped from original
ncid = netcdf.open([spath 'full.ocean_static_FishMIP_remapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% map
clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%%
test = double(squeeze(nmdz(:,:,6)));
test2 = double(squeeze(nlgz(:,:,6)));

geolon = double(geolon);
geolat = double(geolat);
wet = double(wet);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',mlonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,test)
colormap('jet')
colorbar

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',mlonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,wet)
colormap('jet')
colorbar

figure
pcolor(test); shading flat;
colorbar

%%
hpmdz_long_name     = 'medium zooplankton loss to higher preds. integrated in top 100m';
hpmdz_units         = 'gC m-2 d-1';
hplgz_long_name     = 'large zooplankton loss to higher preds.integrated in top 100m';
hplgz_units         = 'gC m-2 d-1';

%%
save([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.mat'],...
    'time','time_units','yr',...
    'hpmdz_long_name','hplgz_long_name','hpmdz_units','hplgz_units',...
    'LAT','LON','hploss_nmdz_100','hploss_nlgz_100','-v7.3');

