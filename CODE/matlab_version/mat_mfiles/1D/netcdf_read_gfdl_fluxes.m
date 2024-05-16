% Read GFDL 1/4 netcdfs
% nitrogen prim. prod. integrals in upper 100m
% Raw units provided by GFDL
% Before conversion by Matthias

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/GFDL_reanalysis/';
fpath='/Volumes/petrik-lab/Fish-MIP/Phase3/GFDL_reanalysis/';

%% one file
ncdisp([fpath '20000101.ocean_cobalt_fluxes_int_FishMIP.nc'])

%%
% Global Attributes:
% filename         = '20000101.ocean_cobalt_fluxes_int.nc'
% title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
% associated_files = 'areacello: 20000101.ocean_static.nc'
% grid_type        = 'regular'
% grid_tile        = 'N/A'
% history          = 'Thu Mar 10 09:35:08 2022: ncks -v jprod_nlgp_100,jprod_ndi_100,jprod_nsmp_100 20000101.ocean_cobalt_fluxes_int.nc 20000101.ocean_cobalt_fluxes_int_FishMIP.nc'
% NCO              = '"4.5.4"'
% missing_value = 1.000000020040877e+20
% _FillValue    = 1.000000020040877e+20

% Dimensions:
% time = 12    (UNLIMITED)
% yh   = 1080
% xh   = 1440

% jprod_ndi_100
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'Diazotroph nitrogen prim. prod. integral in upper 100m'
% units         = 'mol m-2 s-1'

% jprod_nlgp_100
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'Large phyto. nitrogen  prim. prod. integral in upper 100m'
% units         = 'mol m-2 s-1'

% jprod_nsmp_100
% Size:       1440x1080x12
% Dimensions: xh,yh,time
% long_name     = 'Small phyto. nitrogen  prim. prod. integral in upper 100m'
% units         = 'mol m-2 s-1'

% time
% Size:       12x1
% Dimensions: time
% long_name      = 'time'
% units          = 'days since 1959-01-01 00:00:00'
% cartesian_axis = 'T'
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'
% bounds         = 'time_bnds'

% xh
% Size:       1440x1
% Dimensions: xh
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

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time+1)/12)+1601;
runs = find(yr>1950 & yr<=2015);
z100 = find(lev <= 100);

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
    
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_thkcello_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Integrate top 100 m
zmeso_100 = squeeze(nansum((zmeso.*thkcello),3));

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

% save([fpath 'gfdl_hist_zmeso_100_monthly_1950_2014.mat'],'zmeso_100','yr',...
%     'long_name','standard_name','units_orig','units_vint','lat','lon',...
%     'runs','z100','lev');
%




