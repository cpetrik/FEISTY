% Read CMIP6 netcdfs
% IPSL SSP 126
% TB

clear
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';

%% Temp btm
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_tob_60arcmin_global_monthly_2101_2300.nc'])

%%
standard_name      = 'sea_water_potential_temperature';
long_name          = 'Sea Water Potential Temperature';
units              = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_tob_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
tob(tob >= 1.00e+20) = NaN;

% Time
yr = ((time)/12)+1601;

%%
temp_btm = double(tob);
clear tob
ztest2 = squeeze(temp_btm(:,:,10));

figure
pcolor(ztest2); shading flat

%%
save([fpath 'ipsl_ssp126_temp_btm_monthly_2101_2300.mat'],'temp_btm','yr',...
    'long_name','standard_name','units','lat','lon','time');





