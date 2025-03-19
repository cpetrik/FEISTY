% Read CMIP6 netcdfs
% CESM SSP585
% TB

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp585/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp585/';
tpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';

%% zoo zall
ncdisp([fpath 'cesm2-waccm_r1i1p1f1_ssp585_tob_60arcmin_global_monthly_2015_2299.nc'])

%%
standard_name = 'mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Zooplankton Carbon Concentration';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       320x384x15x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1 00:00:00'

%% 
ncid = netcdf.open([fpath 'cesm2-waccm_r1i1p1f1_ssp585_tob_60arcmin_global_monthly_2015_2299.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
tob(tob >= 1.00e+20) = NaN;

%% Get subsets of tob so smaller memory?
yr = ((time+1)/12)+1601;

temp_btm = double(tob);
clear tob 

save([fpath 'cesm2_ssp585_temp_btm_monthly_2015_2299.mat'],'temp_btm',...
    'yr','long_name','standard_name','units','lat','lon',...
    'time');
