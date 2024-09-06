% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss all depths, integrate top 100m
% Regridded by code merging Xiao's with 1 line from Matthias
% Fixed LZ only

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% zmeso regular grid
ncdisp([fpath '19810101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc'])

%Global Attributes:
%            CDI              = 'Climate Data Interface version 2.4.2 (https://mpimet.mpg.de/cdi)'
%            Conventions      = 'CF-1.6'
%            filename         = '19810101.ocean_cobalt_fluxes_int.nc'
%            title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
%            associated_files = 'areacello: 19810101.ocean_static.nc'
%            grid_type        = 'regular'
%            grid_tile        = 'N/A'
%            history          = 'Fri Aug 02 16:48:43 2024: cdo remapbil,r1440x720 -selname,jhploss_nlgz_100 /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19810101.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19810101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc
%                               Fri Aug  2 16:45:37 2024: ncatted -O -a coordinates,jhploss_nlgz_100,c,c,geolat geolon /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19810101.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc
%                               Fri Aug 02 16:44:17 2024: cdo -O -merge /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19810101.ocean_cobalt_fluxes_int_FishMIP_CP.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/full.ocean_static_FishMIP.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19810101.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc
%                               Fri Nov 17 12:22:06 2023: ncks -v jhploss_nmdz_100,jhploss_nlgz_100 19810101.ocean_cobalt_fluxes_int.nc 19810101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'
%            NCO              = 'netCDF Operators version 5.2.7 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco, Citation = 10.1016/j.envsoft.2008.03.004)'
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


%% 
ncid = netcdf.open([fpath '19810101.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
jhploss_nlgz_100(jhploss_nlgz_100 >= 1.00e+19) = NaN;

inan = isnan(jhploss_nlgz_100(:));
ocell = ~isnan(squeeze(jhploss_nlgz_100(:,:,1,1)));
lcell = isnan(squeeze(jhploss_nlgz_100(:,:,1,1)));

%% molN/m2/s --> gC/m2/d
hploss_nlgz_100 = double(jhploss_nlgz_100 * (106/16) * 12.01 * 60*60*24);

%%
test = double(squeeze(hploss_nlgz_100(:,:,6)));

figure
pcolor(test); shading flat; colorbar; 
colormap('jet')

[rLAT,rLON] = meshgrid(lat,lon);

%% Check if same grid as other ISIMIP files =====================
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
load([qpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])
load([qpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])

%% Same orientation?
figure(11)
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

%flipped LR and shifted in matrix
figure(5)
pcolor(fliplr(test)); shading flat; 
colorbar
clim([0 0.01])
title('HP')

figure(6)
pcolor(deptho); shading flat
colorbar
title('depth')

%% Depth has a few more ocean cells
oid = find(~isnan(squeeze(jhploss_nlgz_100(:,:,1,1))));

% oid       670568x1
% GRD.ID    670589x1
gdiff=setdiff(GRD.ID,oid); %off because shifted

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
clim([0 0.05]);

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);



%% Need to shift to use with ISIMIP vars, or match up lat-lon



