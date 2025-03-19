% Read CMIP6 netcdfs
% CESM Hist
% Det

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
tpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';

%% zoo zall
ncdisp([fpath 'cesm2-waccm_r1i1p1f1_historical_expc-bot_60arcmin_global_monthly_1850_2014.nc'])

%%
long_name     = 'Downward Flux of Particulate Organic Carbon on Bottom (z_b)';
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
det_units     = 'mol m-2 s-1';
missing_value = 1.000000020040877e+20;
%Size:       320x384x15x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1 00:00:00'

%% 
ncid = netcdf.open([fpath 'cesm2-waccm_r1i1p1f1_historical_expc-bot_60arcmin_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
expc(expc >= 1.00e+20) = NaN;

%% Get subsets of expc-bot so smaller memory?
yr = ((time+1)/12)+1601;

det = double(expc);
clear expc

save([fpath 'cesm2_hist_det_monthly_1850_2014.mat'],'det',...
    'yr','long_name','standard_name','det_units','lat','lon',...
    'time');
