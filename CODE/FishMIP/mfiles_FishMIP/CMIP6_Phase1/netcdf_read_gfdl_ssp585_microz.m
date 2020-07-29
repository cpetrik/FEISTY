% Read Fish-MIP netcdfs
% GFDL ssp585

clear all
close all

fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/ssp585/';

%% zmicro
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp585_zmicro_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_zmicroorophyll_in_sea_water';
long_name     = 'Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water';
units         = 'kg m-3';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp585_zmicro_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

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
%zmicro: 360x180x35x1032
%time: 1980
%NaNs = 1.0000e+20

%% Get subset of zmicro
% Time 
yr = ((time)/12)+1601-1;
z100 = find(lev <= 100);

for i = nvars
    varname = netcdf.inqVar(ncid, i-1);
    zmicro = netcdf.getVar(ncid,i-1, [0,0,0,0],[360 180 length(z100) length(time)]);
    zmicro(zmicro >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m & subset time
zmicro_100 = squeeze(nansum(zmicro,3));

%% 
save([fpath 'gfdl_ssp585_zmicro100_monthly_2015_2100.mat'],'zmicro_100','time',...
    'yr','long_name','standard_name','units','lev','z100');





