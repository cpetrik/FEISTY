% Read GFDL 1/4 netcdfs
% Horiz Velocities and Salinity 3D
% Raw units provided by GFDL
% Before conversion by Matthias

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_month_z_FishMIP.nc'])

%%
% Global Attributes:
% filename         = '20000101.ocean_month_z.nc'
% title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% associated_files = 'areacello: 20000101.ocean_static.nc'
% grid_type        = 'regular'
% grid_tile        = 'N/A'
% history          = 'Wed Mar  9 22:38:21 2022: ncks -v uo,vo,so 20000101.ocean_month_z.nc 20000101.ocean_month_z_FishMIP.nc'
% NCO              = '"4.5.4"'

% Dimensions:
% time = 12    (UNLIMITED)
% z_l  = 35
% yh   = 1080
% xh   = 1440
% xq   = 1440
% yq   = 1080

% so
% Size:       1440x1080x35x12
% Dimensions: xh,yh,z_l,time
% Datatype:   single
% Attributes:
% long_name     = 'Sea Water Salinity'
% units         = 'psu'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'
% standard_name = 'sea_water_salinity'

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

% uo
% Size:       1440x1080x35x12
% Dimensions: xq,yh,z_l,time
% Datatype:   single
% Attributes:
% long_name     = 'Sea Water X Velocity'
% units         = 'm s-1'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'z_l:mean yh:mean xq:point time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% standard_name = 'sea_water_x_velocity'
% interp_method = 'none'

% vo
% Size:       1440x1080x35x12
% Dimensions: xh,yq,z_l,time
% Datatype:   single
% Attributes:
% long_name     = 'Sea Water Y Velocity'
% units         = 'm s-1'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'z_l:mean yq:point xh:mean time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% standard_name = 'sea_water_y_velocity'
% interp_method = 'none'

% xh
% Size:       1440x1
% Dimensions: xh
% Datatype:   double
% Attributes:
% long_name      = 'h point nominal longitude'
% units          = 'degrees_east'
% cartesian_axis = 'X'

% xq
% Size:       1440x1
% Dimensions: xq
% Datatype:   double
% Attributes:
% long_name      = 'q point nominal longitude'
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

% yq
% Size:       1080x1
% Dimensions: yq
% Datatype:   double
% Attributes:
% long_name      = 'q point nominal latitude'
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




