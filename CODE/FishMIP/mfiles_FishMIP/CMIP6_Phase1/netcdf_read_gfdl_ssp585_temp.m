% Read Fish-MIP netcdfs
% GFDL ssp585

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/ssp585/';
fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/ssp585/';

%% thetao
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp585_thetao_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp585_thetao_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%lev: 35
%thetao: 360x180x35x1032
%time: 1032
%NaNs = 1.0000e+20

%% Integrate top 100 m 
z100 = find(lev <= 100);
temp_100 = squeeze(nanmean(thetao(:,:,z100,:),3));

%% Time 
yr = ((time)/12)+1601-1;

save([fpath 'gfdl_ssp585_temp100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','lev','z100');
save([spath 'gfdl_ssp585_temp100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','lev','z100');





