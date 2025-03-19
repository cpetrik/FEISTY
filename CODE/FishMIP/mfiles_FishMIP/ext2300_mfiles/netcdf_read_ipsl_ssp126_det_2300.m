% Read CMIP6 netcdfs
% IPSL SSP 126
% Det

clear
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';

%% Det btm
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_expc-bot_60arcmin_global_monthly_2101_2300.nc'])

%%
long_name     = 'Sinking Particulate Organic Carbon Flux on Bottom (z_b)';
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
units         = 'mol m-2 s-1';
missing_value = 1.000000020040877e+20;
%Size:       360x180x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_expc-bot_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
expc(expc >= 1.00e+20) = NaN;

% Time
yr = ((time)/12)+1601;

%%
det = double(expc);
clear expc
ztest2 = squeeze(det(:,:,10))*12.01*9*60*60*24;

figure
pcolor(ztest2); shading flat
colorbar
clim([0 3])

%%
save([fpath 'ipsl_ssp126_det_monthly_2101_2300.mat'],'det','yr',...
    'long_name','standard_name','units','lat','lon',...
    'time');





