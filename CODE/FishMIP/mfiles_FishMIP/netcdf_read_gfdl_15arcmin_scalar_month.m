% Read GFDL 1/4 netcdfs
% 12 month vector
% Net Downward Shortwave Radiation at Sea Water Surface Area Averaged
% Raw units provided by GFDL
% Before conversion by Matthias

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_scalar_month_FishMIP.nc'])

%%
% Global Attributes:
% filename  = '20000101.ocean_scalar_month.nc'
% title     = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% grid_type = 'regular'
% grid_tile = 'N/A'
% history   = 'Thu Mar 10 09:30:33 2022: ncks -v ave_rsntds 20000101.ocean_scalar_month.nc 20000101.ocean_scalar_month_FishMIP.nc'
% NCO       = '"4.5.4"'

% Dimensions:
% time        = 12    (UNLIMITED)
% scalar_axis = 1
% Variables:

% ave_rsntds
% Size:       1x12
% Dimensions: scalar_axis,time
% Datatype:   single
% Attributes:
% long_name     = 'Net Downward Shortwave Radiation at Sea Water Surface Area Averaged'
% units         = 'W m-2'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20
% cell_methods  = 'time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% standard_name = 'net_downward_shortwave_flux_at_sea_water_surface_area_averaged'

% scalar_axis
% Size:       1x1
% Dimensions: scalar_axis
% Datatype:   double
% Attributes:
% long_name      = 'none'
% units          = 'none'
% cartesian_axis = 'N'

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




