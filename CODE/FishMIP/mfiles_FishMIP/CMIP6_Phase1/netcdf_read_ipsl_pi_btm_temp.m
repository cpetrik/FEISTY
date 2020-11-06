% Read Fish-MIP netcdfs
% bottom temp

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/preindust/';

%%
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_tob_onedeg_global_monthly_1601_2100.nc'])

standard_name = 'sea_water_potential_temperature_at_sea_floor';
long_name     = 'Sea Water Potential Temperature at Sea Floor';
units         = 'degC';
%comment       = 'Model data on the 1x1 grid includes values in all cells for which any ocean exists on the native grid. For mapping purposes, we recommend using a land mask such as World Ocean Atlas to cover these areas of partial land.  For calculating approximate integrals, we recommend multiplying by cell area (areacello).'

%time = months since 1601-1-1 00:00:00

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_tob_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%NaNs on land cells
tob(tob >= 1.0000e+20) = NaN;

%% Different time periods
yr = ((time+1)/12)+1601; %months since 1601-1-1 = first month is Jan 1601

spin = find(yr>1850 & yr<=1950);
runs = find(yr>1950);

%% spinup
temp_btm = tob(:,:,spin);
y = yr(spin);
save([fpath 'ipsl_pi_temp_btm_monthly_1850_1949.mat'],'temp_btm','time',...
    'yr','spin','y','long_name','standard_name','units');

%% regular runs
clear temp_btm y
temp_btm = tob(:,:,runs);
y = yr(runs);
save([fpath 'ipsl_pi_temp_btm_monthly_1950_2100.mat'],'temp_btm','time',...
    'yr','runs','y','long_name','standard_name','units');




