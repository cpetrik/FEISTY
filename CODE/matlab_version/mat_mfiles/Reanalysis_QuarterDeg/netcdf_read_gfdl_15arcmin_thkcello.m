% Read GFDL 1/4 netcdfs
% ctrl thkcello

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath 'gfdl-mom6-cobalt2_ctrlclim_thkcello_15arcmin_global_fixed.nc'])

%%
% Source:
%            /Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/gfdl-mom6-cobalt2_ctrlclim_thkcello_15arcmin_global_fixed.nc
% Format:
%            netcdf4_classic
% Global Attributes:
%            CDI                      = 'Climate Data Interface version 1.9.8 (https://mpimet.mpg.de/cdi)'
%            Conventions              = 'CF-1.6'
%            cdo_openmp_thread_number = 4
%            isimip_id                = 'd51332d6-b45c-40be-871f-11e4fdd4fddc'
% Dimensions:
%            time = 1     (UNLIMITED)
%            lon  = 1440
%            lat  = 720
%            lev  = 35
% Variables:
%     time    
%            Size:       1x1
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
%     lev     
%            Size:       35x1
%            Dimensions: lev
%            Datatype:   double
%            Attributes:
%                        long_name      = 'Depth at cell center'
%                        units          = 'meters'
%                        positive       = 'down'
%                        axis           = 'Z'
%                        cartesian_axis = 'Z'
%                        edges          = 'z_i'
%     thkcello
%            Size:       1440x720x35x1
%            Dimensions: lon,lat,lev,time
%            Datatype:   single
%            Attributes:
%                        standard_name = 'cell_thickness'
%                        long_name     = 'Sea Water Potential Temperature'
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20

%%
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_ctrlclim_thkcello_15arcmin_global_fixed.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end





