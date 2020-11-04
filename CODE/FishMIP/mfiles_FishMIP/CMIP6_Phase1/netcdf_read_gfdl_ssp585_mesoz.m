% Read Fish-MIP netcdfs
% GFDL ssp585

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/ssp585/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';

%% zmeso
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp585_zmeso_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton expressed as Carbon in sea water';
units         = 'mol m-3';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp585_zmeso_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Vars
%lat: 180
%lon: 360
%lev: 35
%zmeso: 360x180x35x1032
%time: 1980
%NaNs = 1.0000e+20

%% Get subset of zmeso
% Time 
yr = ((time)/12)+1601-1;
z100 = find(lev <= 100);

for i = nvars
    varname = netcdf.inqVar(ncid, i-1);
    zmeso = netcdf.getVar(ncid,i-1, [0,0,0,0],[360 180 length(z100) length(time)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m & subset time
zmeso_100 = squeeze(nansum(zmeso,3));

%% 
save([fpath 'gfdl_ssp585_zmeso100_monthly_2015_2100.mat'],'zmeso_100','time',...
    'yr','long_name','standard_name','units','lev','z100');
save([spath 'gfdl_ssp585_zmeso100_monthly_2015_2100.mat'],'zmeso_100','time',...
    'yr','long_name','standard_name','units','lev','z100');





