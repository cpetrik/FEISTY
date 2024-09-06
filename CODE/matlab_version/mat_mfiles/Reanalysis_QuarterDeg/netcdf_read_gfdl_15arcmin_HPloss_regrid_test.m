% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss all depths, integrate top 100m
% 1/4 grid from ISIMIP

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% zmeso native grid
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


%% zmeso regular grid
ncdisp([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc'])

% Source:
%            /Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/19610101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc
% Format:
%            64bit
% Global Attributes:
%            CDI              = 'Climate Data Interface version 2.4.2 (https://mpimet.mpg.de/cdi)'
%            Conventions      = 'CF-1.6'
%            filename         = '19610101.ocean_cobalt_fluxes_int.nc'
%            title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
%            associated_files = 'areacello: 19610101.ocean_static.nc'
%            grid_type        = 'regular'
%            grid_tile        = 'N/A'
%            history          = 'Thu Jul 18 11:57:03 2024: cdo -selname,jhploss_nlgz_100,jhploss_nmdz_100 -sellonlatbox,-180,180,-90,90 -remapbil,grid.des 19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc.tmp 19610101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc
%                               Thu Jul 18 11:57:02 2024: cdo -O -merge 19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc full.ocean_static_FishMIP.nc 19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc.tmp
%                               Fri Nov 17 12:22:25 2023: ncks -v jhploss_nmdz_100,jhploss_nlgz_100 19610101.ocean_cobalt_fluxes_int.nc 19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'
%            NCO              = 'netCDF Operators version 5.0.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)'
%            CDO              = 'Climate Data Operators version 2.4.2 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            time = 12    (UNLIMITED)
%            lon  = 1440
%            lat  = 720
% Variables:
%     time            
%            Size:       12x1
%            Dimensions: time
%            Datatype:   double
%            Attributes:
%                        standard_name = 'time'
%                        long_name     = 'time'
%                        units         = 'days since 1959-01-01 00:00:00'
%                        calendar      = 'standard'
%                        axis          = 'T'
%     lon             
%            Size:       1440x1
%            Dimensions: lon
%            Datatype:   double
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude'
%                        units         = 'degrees_east'
%                        axis          = 'X'
%     lat             
%            Size:       720x1
%            Dimensions: lat
%            Datatype:   double
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude'
%                        units         = 'degrees_north'
%                        axis          = 'Y'
%     jhploss_nlgz_100
%            Size:       1440x720x12
%            Dimensions: lon,lat,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m'
%                        units         = 'mol m-2 s-1'
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'
%     jhploss_nmdz_100
%            Size:       1440x720x12
%            Dimensions: lon,lat,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m'
%                        units         = 'mol m-2 s-1'
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'



%% 
ncid = netcdf.open([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc'],'NC_NOWRITE');
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
jhploss_nmdz_100(jhploss_nmdz_100 >= 1.00e+19) = NaN;
jhploss_nlgz_100(jhploss_nlgz_100 >= 1.00e+19) = NaN;

inan = isnan(jhploss_nmdz_100(:));
ocell = ~isnan(squeeze(jhploss_nmdz_100(:,:,1,1)));
lcell = isnan(squeeze(jhploss_nmdz_100(:,:,1,1)));

%% molN/m2/s --> gC/m2/d
hploss_nmdz_100 = double(jhploss_nmdz_100 * (106/16) * 12.01 * 60*60*24);
hploss_nlgz_100 = double(jhploss_nlgz_100 * (106/16) * 12.01 * 60*60*24);

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

figure
pcolor(geolon_c,geolat_c,test); shading flat; colorbar; 
colormap('jet')

figure
pcolor(geolon_c,geolat_c,test2); shading flat; colorbar;
colormap('jet')

geolon_c = double(geolon_c);
geolat_c = double(geolat_c);
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',mlonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_c,geolon_c,test)
colormap('jet')
colorbar


%% Check if same grid as other ISIMIP files =====================
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
load([qpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])
load([qpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])

%% Same orientation?
figure(1)
pcolor(rLAT); shading flat;
colorbar
title('rLAT')

figure(2)
pcolor(LAT); shading flat;
colorbar
title('LAT')

figure(3)
pcolor(rLON); shading flat;
colorbar
title('LON')

figure(4)
pcolor(LON); shading flat;
colorbar
title('LON')

figure(5)
pcolor(test); shading flat;
colorbar
clim([0 0.1])
title('HP')

figure(6)
pcolor(deptho); shading flat
colorbar
title('depth')

%% Depth has more ocean cells
oid = find(~isnan(squeeze(jhploss_nmdz_100(:,:,1,1))));

% oid       640165x1
% GRD.ID    670589x1
gdiff=setdiff(GRD.ID,oid); %40918

[ni,nj] = size(LON);
gMAT = nan*ones(ni,nj);
gMAT(gdiff) = 1;

figure(7)
pcolor(gMAT); shading flat
title('ocean cell is depth and not remapped')

%% map
latlim=[-90 90];
lonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(rLAT,rLON,test)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);


figure
axesm('stereo','MapLatLimit',[65 90])
surfm(rLAT,rLON,test)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);

figure
axesm('stereo','MapLatLimit',[65 90])
surfm(LAT,LON,deptho)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(LAT,LON,gMAT)

figure
axesm('stereo','MapLatLimit',[65 90])
surfm(LAT,LON,gMAT)

%% What about combined mesozoo orig, missing data?
gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

ncdisp([gpath '20000101.ocean_cobalt_fluxes_int_FishMIP.nc'])

%% 
ncid = netcdf.open([gpath '20000101.ocean_cobalt_fluxes_int_FishMIP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

[oLAT,oLON] = meshgrid(yh,xh);

%%
jprod_diat_100 = double(jprod_ndi_100 * (106/16) * 12.01 * 60*60*24);
test3 = (squeeze(jprod_diat_100(:,:,6)));


figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(oLAT,oLON,test3)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.05]);

figure
axesm('stereo','MapLatLimit',[65 90])
surfm(oLAT,oLON,test3)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.05]);











