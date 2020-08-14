% Read Fish-MIP netcdfs
% bottom temp

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/IPSL/hist/';

%%
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_tob_onedeg_global_monthly_1850_2014.nc'])

standard_name = 'sea_water_potential_temperature_at_sea_floor';
long_name     = 'Sea Water Potential Temperature at Sea Floor';
units         = 'degC';
%comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'

%time = months since 1601-1-1 00:00:00

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_tob_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
tob(tob >= 1.0000e+20) = NaN;

%% Time
%yr = ((time+1)/12)+1601-1;
yr = ((time)/12)+1601-1;
runs = find(yr>1949 & yr<=2015); 
temp_btm = tob(:,:,runs);

save([fpath 'ipsl_hist_temp_btm_monthly_1950_2014.mat'],'temp_btm','time',...
    'yr','runs','long_name','standard_name','units');



