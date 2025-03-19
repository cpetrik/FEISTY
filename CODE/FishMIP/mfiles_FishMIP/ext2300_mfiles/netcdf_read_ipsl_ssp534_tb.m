% Read CMIP6 netcdfs
% IPSL SSP 534
% Temp Btm

clear
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';

%% Temp zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp534-over_thetao-bot_60arcmin_global_monthly_2040_2300.nc'])

%%
long_name     = 'Sea Water Potential Temperature on Bottom (z_b)';
standard_name = 'sea_water_potential_temperature';
units         = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp534-over_thetao-bot_60arcmin_global_monthly_2040_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time+1)/12)+1601;

thetao(thetao >= 1.00e+20) = NaN;

%%
temp_btm = double(thetao);
clear thetao
ztest2 = squeeze(temp_btm(:,:,100));

figure
pcolor(ztest2)
shading flat

%%
save([fpath 'ipsl_ssp534-over_temp_btm_monthly_2040_2300.mat'],'temp_btm','yr',...
    'long_name','standard_name','units','lat','lon','time');





