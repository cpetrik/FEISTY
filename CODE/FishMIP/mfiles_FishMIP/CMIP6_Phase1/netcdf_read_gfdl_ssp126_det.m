% Read Fish-MIP netcdfs
% bottom POC flux

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/ssp126/';

%%
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp126_expc-bot_onedeg_global_monthly_2015_2100.nc'])

long_name     = 'Sinking Particulate Organic Carbon Flux on Bottom (z_b)';
standard_name = 'sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water';
units         = 'mol m-2 s-1';
%comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'

%time = months since 1601-1-1 00:00:00

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp126_expc-bot_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
expc(expc >= 1.0000e+20) = NaN;

%% Time 
yr = ((time)/12)+1601-1;
det_btm = expc;

save([fpath 'gfdl_ssp126_det_btm_monthly_2015_2100.mat'],'det_btm','time',...
    'yr','long_name','standard_name','units');




