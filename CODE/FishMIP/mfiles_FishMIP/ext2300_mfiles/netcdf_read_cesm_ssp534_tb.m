% Read CMIP6 netcdfs
% CESM SSP534-over
% TB

clear
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';
tpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';

%% zoo zall
ncdisp([fpath 'cesm2-waccm_r1i1p1f1_ssp534-over_thetao-bot_60arcmin_global_monthly_2040_2299.nc'])

%%
long_name     = 'Sea Water Potential Temperature on Bottom (z_b)';
standard_name = 'sea_water_potential_temperature';
units         = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       320x384x15x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1 00:00:00'

%%
ncid = netcdf.open([fpath 'cesm2-waccm_r1i1p1f1_ssp534-over_thetao-bot_60arcmin_global_monthly_2040_2299.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
thetao(thetao >= 1.00e+20) = NaN;

%% 
yr = ((time+1)/12)+1601;

temp_btm = double(thetao);
clear thetao

save([fpath 'cesm2_ssp534-over_temp_btm_monthly_2040_2299.mat'],'temp_btm',...
    'yr','long_name','standard_name','units','lat','lon',...
    'time');
