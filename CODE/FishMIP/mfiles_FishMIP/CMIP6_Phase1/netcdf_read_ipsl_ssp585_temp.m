% Read Fish-MIP netcdfs
% IPSL ssp585

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/ssp585/';
fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/IPSL/ssp585/';

%% thetao
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_thetao_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_thetao_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Vars
%lat: 180
%lon: 360
%lev: 75
%thetao: 360x180x75x1032
%time: 1032
%NaNs = 1.0000e+20

%% Get subset of thetao
% Time
yr = ((time)/12)+1601-1;
runs = find(yr>2014 & yr<=2100);
z100 = find(olevel <= 100);

for i = nvars
    varname = netcdf.inqVar(ncid, i-1);
    thetao = netcdf.getVar(ncid,i-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
    thetao(thetao >= 1.00e+20) = NaN;
end
netcdf.close(ncid);

%% Integrate top 100 m 
temp_100 = squeeze(nanmean(thetao,3));

%% Time 

save([fpath 'ipsl_ssp585_temp100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','olevel','z100');
save([spath 'ipsl_ssp585_temp100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','olevel','z100');





