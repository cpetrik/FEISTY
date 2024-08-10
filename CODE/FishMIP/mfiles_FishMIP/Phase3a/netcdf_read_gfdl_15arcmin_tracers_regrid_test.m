% Read GFDL netcdfs
% obsclim
% mdz & lgz biomass all depths
% Regridded by code merging Xiao's with 1 line from Matthias
% works on both kelpfish and tunafish now
% needed to add nco tools to kelpfish

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% zmeso regular grid
ncdisp([fpath '19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc'])

% Global Attributes:
% CDI              = 'Climate Data Interface version 2.4.2 (https://mpimet.mpg.de/cdi)'
%            Conventions      = 'CF-1.6'
%            filename         = '19710101.ocean_cobalt_tracers_month_z.nc'
%            title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
%            associated_files = 'areacello: 19710101.ocean_static.nc'
%            grid_type        = 'regular'
%            grid_tile        = 'N/A'
%            history          = 'Thu Aug 08 12:05:43 2024: cdo remapbil,r1440x720 -selname,nlgz,nmdz /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_temp.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc
%                               Thu Aug  8 10:40:44 2024: ncatted -O -a coordinates,nmdz,c,c,geolat geolon /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_temp.nc
%                               Thu Aug  8 10:15:22 2024: ncatted -O -a coordinates,nlgz,c,c,geolat geolon /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_temp.nc
%                               Thu Aug 08 09:59:21 2024: cdo -O -merge /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19710101.ocean_cobalt_tracers_month_z_FishMIP_CP.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/full.ocean_static_FishMIP.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_temp.nc
%                               Fri Nov 17 12:19:14 2023: ncks -v nmdz,nlgz 19710101.ocean_cobalt_tracers_month_z.nc 19710101.ocean_cobalt_tracers_month_z_FishMIP_CP.nc'
%            NCO              = 'netCDF Operators version 5.2.7 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco, Citation = 10.1016/j.envsoft.2008.03.004)'
%            CDO              = 'Climate Data Operators version 2.4.2 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            time = 12    (UNLIMITED)
%            lon  = 1440
%            lat  = 720
%            z_l  = 35
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
%     z_l 
%            Size:       35x1
%            Dimensions: z_l
%            Datatype:   double
%            Attributes:
%                        long_name      = 'Depth at cell center'
%                        units          = 'meters'
%                        positive       = 'down'
%                        axis           = 'Z'
%                        cartesian_axis = 'Z'
%                        edges          = 'z_i'
%     nlgz
%            Size:       1440x720x35x12
%            Dimensions: lon,lat,z_l,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'large Zooplankton Nitrogen'
%                        units         = 'mol/kg'
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'
%     nmdz
%            Size:       1440x720x35x12
%            Dimensions: lon,lat,z_l,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Medium-sized zooplankton Nitrogen'
%                        units         = 'mol/kg'
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'

%% 
ncid = netcdf.open([fpath '19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-2)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%%
for i = (nvars-1):nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
nmdz(nmdz >= 1.00e+19) = NaN;
nlgz(nlgz >= 1.00e+19) = NaN;

%%
testM = (squeeze(nmdz(:,:,1,6)));
testL = (squeeze(nlgz(:,:,1,7)));

figure
pcolor(testM); shading flat; colorbar; 

figure
pcolor(testL); shading flat; colorbar; 

%% Check if same grid as other ISIMIP files =====================
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
load([qpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])

%% Need to shift to use with ISIMIP vars, or match up lat-lon
%Xiao: split, flip left and right, and concatenate the matrix before plotting the map

Mtsplit1 = testM(1:720,:);
Mtsplit2 = testM(721:end,:);
Mtflip1 = fliplr(Mtsplit1);
Mtflip2 = fliplr(Mtsplit2);
Mtcomb = [Mtflip2;Mtflip1]; %This looks aligned with ISIMIP depth

Ltsplit1 = testL(1:720,:);
Ltsplit2 = testL(721:end,:);
Ltflip1 = fliplr(Ltsplit1);
Ltflip2 = fliplr(Ltsplit2);
Ltcomb = [Ltflip2;Ltflip1]; %This looks aligned with ISIMIP depth

close all

figure
pcolor(Mtcomb); shading flat

figure
pcolor(Ltcomb); shading flat

figure
pcolor(deptho); shading flat




