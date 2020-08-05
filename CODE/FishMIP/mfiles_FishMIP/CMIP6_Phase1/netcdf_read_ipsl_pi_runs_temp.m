% Read Fish-MIP netcdfs
% IPSL PreIndust

clear all
close all

fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/preindust/';
spath='/Volumes/FEISTY/Fish-MIP/CMIP6/IPSL/preindust/';

%% Temp zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_thetao_onedeg_global_monthly_1601_2100.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_thetao_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');

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
%thetao: 360x180x75x6000
%time: 6000
%NaNs = 1.0000e+20

%% Get subset of temp
% Time
yr = ((time+1)/12)+1601-1;
runs = find(yr>1950 & yr<=2100);
z100 = find(olevel <= 100);

i = nvars;
thetao = netcdf.getVar(ncid,i-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
thetao(thetao >= 1.00e+20) = NaN;
netcdf.close(ncid);

%% Mean top 100 m 
temp_100 = squeeze(nanmean(thetao,3));

%%
save([fpath 'ipsl_pi_temp100_monthly_1950_2100.mat'],'temp_100','time',...
    'yr','runs','long_name','standard_name','units','olevel','z100');
save([spath 'ipsl_pi_temp100_monthly_1950_2100.mat'],'temp_100','time',...
    'yr','runs','long_name','standard_name','units','olevel','z100');





