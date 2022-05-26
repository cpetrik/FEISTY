% Read GFDL 1/4 netcdfs
% Thetao 3D
% Raw units provided by GFDL
% Before conversion by Matthias

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_month_z_FishMIP_b.nc'])

%%
% Global Attributes:
% filename         = '20000101.ocean_month_z.nc'
% title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% associated_files = 'areacello: 20000101.ocean_static.nc'
% grid_type        = 'regular'
% grid_tile        = 'N/A'
% history          = 'Mon Mar 21 09:36:08 2022: ncks -v thetao 20000101.ocean_month_z.nc 20000101.ocean_month_z_FishMIP_b.nc'
% NCO              = '"4.5.4"'
% Dimensions:
% time = 12    (UNLIMITED)
% z_l  = 35
% yh   = 1080
% xh   = 1440

% thetao
% Size:       1440x1080x35x12
% Dimensions: xh,yh,z_l,time
% Datatype:   single
% Attributes:
% long_name     = 'Sea Water Potential Temperature'
% units         = 'degC'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'
% standard_name = 'sea_water_potential_temperature'

% time
% Size:       12x1
% Dimensions: time
% Datatype:   double
% Attributes:
% long_name      = 'time'
% units          = 'days since 1959-01-01 00:00:00'
% cartesian_axis = 'T'
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'
% bounds         = 'time_bnds'

% xh
% Size:       1440x1
% Dimensions: xh
% Datatype:   double
% Attributes:
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
% cartesian_axis = 'Z'
% positive       = 'down'
% edges          = 'z_i'

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time




% save([fpath 'gfdl_hist_zmeso_100_monthly_1950_2014.mat'],'zmeso_100','yr',...
%     'long_name','standard_name','units_orig','units_vint','lat','lon',...
%     'runs','z100','lev');
%




