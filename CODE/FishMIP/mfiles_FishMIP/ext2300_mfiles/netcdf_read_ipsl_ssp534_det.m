% Read CMIP6 netcdfs
% IPSL SSP 534
% Vert mean top 200 m

clear
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';

%% Det btm
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp534-over_expc-bot_60arcmin_global_monthly_2040_2300.nc'])

%%
long_name     = 'Sinking Particulate Organic Carbon Flux on Bottom (z_b)';
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
units         = 'mol m-2 s-1';
missing_value = 1.000000020040877e+20;
%Size:       360x180x3132
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp534-over_expc-bot_60arcmin_global_monthly_2040_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time+1)/12)+1601;
% runs = find(yr <=2100);
% z200 = find(olevel <= 200);

expc(expc >= 1.00e+20) = NaN;

%%
ztest2 = squeeze(expc(:,:,100)) *12.01*9*60*60*24;

figure
pcolor(ztest2)
shading flat

%%
save([fpath 'ipsl_ssp534-over_det_monthly_2040_2300.mat'],'expc','yr',...
    'long_name','standard_name','units','lat','lon','time');





