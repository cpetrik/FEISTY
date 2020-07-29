% Read Fish-MIP netcdfs
% bottom temp

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/ssp585/';

%%
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp585_thetao-bot_onedeg_global_monthly_2015_2100.nc'])

long_name     = 'Sea Water Potential Temperature on Bottom (z_b)';
standard_name = 'sea_water_potential_temperature';
units         = 'degC';
%comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'

%time = months since 1601-1-1 00:00:00

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp585_thetao-bot_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
thetao(thetao >= 1.0000e+20) = NaN;

%% Time
yr = ((time)/12)+1601-1;
temp_btm = thetao;

save([fpath 'gfdl_ssp585_temp_btm_monthly_2015_2100.mat'],'temp_btm','time',...
    'yr','long_name','standard_name','units');




