% Read Fish-MIP netcdfs
% GFDL Hist

clear all
close all

fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/hist/';

%% chl
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_chl_onedeg_global_monthly_1850_2014.nc'])
standard_name = 'mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water';
long_name     = 'Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water';
units         = 'kg m-3';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_chl_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

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
%chl: 360x180x35x1980
%time: 1980
%NaNs = 1.0000e+20

%% Get subset of chl
% Time 
yr = ((time)/12)+1601-1;
runs = find(yr>1950 & yr<=2014);
z100 = find(lev <= 100);

for i = nvars
    varname = netcdf.inqVar(ncid, i-1);
    chl = netcdf.getVar(ncid,i-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
    chl(chl >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m & subset time
%chl_100 = squeeze(sum(chl(:,:,z100,:),3));
chl_100 = squeeze(sum(chl,3));

%%
save([fpath 'gfdl_hist_chl100_monthly_1950_2014.mat'],'chl_100','time',...
    'yr','runs','long_name','standard_name','units','lev','z100');





