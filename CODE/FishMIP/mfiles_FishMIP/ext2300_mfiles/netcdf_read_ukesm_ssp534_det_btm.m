% Read CMIP6 netcdfs
% UKESM SSP 534
% Det btm

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';

%% det btm
ncdisp([fpath 'ukesm1-0-ll_r4i1p1f2_ssp534-over_expc-bot_60arcmin_global_monthly_2040_2100.nc'])

%%
long_name     = 'Downward Flux of Particulate Organic Carbon on Bottom (z_b)';
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
units         = 'mol m-2 s-1';
missing_value = 1.000000020040877e+20;
%Size:       360x180x732
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
%% 
ncid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp534-over_expc-bot_60arcmin_global_monthly_2040_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
expc(expc >= 1.00e+20) = NaN;

%% Time
yr = ((time)/12)+1601;

ztest = squeeze(expc(:,:,17))*12.01*9*60*60*24;

figure
pcolor(ztest); shading flat
colorbar
clim([0 3])

%%
det = double(expc);
clear expc 

save([fpath 'ukesm_ssp534_det_monthly_2040_2100.mat'],'det','yr',...
    'long_name','standard_name','units','lat','lon','time');





