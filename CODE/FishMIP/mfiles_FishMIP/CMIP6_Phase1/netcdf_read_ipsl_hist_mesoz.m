% Read Fish-MIP netcdfs
% IPSL Hist 1950-2014

clear all
close all

fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/hist/';
spath='/Volumes/FEISTY/Fish-MIP/CMIP6/IPSL/hist/';

%% Meso Zoop zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton expressed as Carbon in sea water';
units         = 'mol m-3';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_zmeso_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

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
%olevel: 75
%zmeso: 360x180x75x1980
%time: 1980
%NaNs = 1.0000e+20

%% Get subset of zmeso
% Time
%yr = ((time+1)/12)+1601-1;
yr = ((time)/12)+1601-1;
runs = find(yr>1949 & yr<=2015);
z100 = find(olevel <= 100);

for i = nvars
    varname = netcdf.inqVar(ncid, i-1);
    zmeso = netcdf.getVar(ncid,i-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
    zmeso(zmeso >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m & subset time
zmeso_100 = squeeze(nansum(zmeso,3));

%%
save([fpath 'ipsl_hist_zmeso100_monthly_1950_2014.mat'],'zmeso_100','time',...
    'yr','runs','long_name','standard_name','units','olevel','z100');
save([spath 'ipsl_hist_zmeso100_monthly_1950_2014.mat'],'zmeso_100','time',...
    'yr','runs','long_name','standard_name','units','olevel','z100');





