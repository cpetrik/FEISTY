% Read Fish-MIP netcdfs
% GFDL ssp126

clear all
close all

fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/ssp126/';
spath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/ssp126/';

%% thetao
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp126_thetao_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp126_thetao_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
%netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%lev: 35
%thetao: 360x180x35x1032
%time: 1032
%NaNs = 1.0000e+20

%% Get subset of temp
% Time 
yr = ((time)/12)+1601-1;
z100 = find(lev <= 100);

for i = nvars
    varname = netcdf.inqVar(ncid, i-1);
    thetao = netcdf.getVar(ncid,i-1, [0,0,0,0],[360 180 length(z100) length(time)]);
    thetao(thetao >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m 
% z100 = find(lev <= 100);
% temp_100 = squeeze(nanmean(thetao(:,:,z100,:),3));

temp_100 = squeeze(nanmean(thetao,3));

%% Time 
yr = ((time)/12)+1601-1;

save([fpath 'gfdl_ssp126_temp100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','lev','z100');
save([spath 'gfdl_ssp126_temp100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','lev','z100');





