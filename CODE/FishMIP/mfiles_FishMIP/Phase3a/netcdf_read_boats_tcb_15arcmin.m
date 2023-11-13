% Read BOATS tcb
% forced with GFDL 1/4 netcdfs
% obsclim

clear 
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/BOATS/';

%% one file
ncdisp([fpath 'boats_gfdl-mom6-cobalt2_obsclim_histsoc_default_tcb_global_monthly_1961_2010.nc'])

%%
% Format: netcdf4_classic

% Global Attributes:
% contact                 = 'Daniele Bianchi <dbianchi@atmos.ucla.edu>, Eric Galbraith <eric.galbraith@mcgill.ca>, Jerome Guiet <jerome.c.guiet@gmail.com>'
% institution             = 'AOS - University of California Los Angeles'
% comment                 = 'Impact model output for ISIMIP3a'
% date_created            = '15/02/23-15:27'
% isimip_qc_pass_date     = '2023-05-22 11:07 UTC'
% isimip_qc_version       = '3.0.2'

% Dimensions:
% lon  = 1440
% lat  = 720
% time = 600   (UNLIMITED)

% Variables:
time_units = 'days since 1901-01-01 00:00:00';
% calendar      = 'gregorian'
% axis          = 'T'
% tcb
% Size:       1440x720x600
% Dimensions: lon,lat,time
% Datatype:   single
% Attributes:
% _FillValue    = 1.000000020040877e+20
% long_name     = 'Total Consumer Biomass Density'
tcb_units = 'g m-2';
% missing_value = 1e+20


%%
ncid = netcdf.open([fpath 'boats_gfdl-mom6-cobalt2_obsclim_histsoc_default_tcb_global_monthly_1961_2010.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

%%
tcb(tcb>1e19) = nan;
tcb = double(tcb);

%% Time
yr = 1901 + (time/365);

%%
save([fpath 'QuarterDeg/gfdl-mom6-cobalt2_ctrlclim_expc-bot_15arcmin_global_monthly_1961_2010.mat'],...
    'GFDL_variable','long_name','standard_name','missing_value','units',...
    'lat','lon','time','LAT','LON','tcb','yr', '-v7.3');

