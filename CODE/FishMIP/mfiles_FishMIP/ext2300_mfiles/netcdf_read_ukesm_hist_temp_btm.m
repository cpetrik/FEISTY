% Read CMIP6 netcdfs
% UKESM Hist
% TB

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

%% temp btm
ncdisp([fpath 'ukesm1-0-ll_r1i1p1f2_historical_tob_60arcmin_global_monthly_1850_2014.nc'])

%%
standard_name      = 'sea_water_potential_temperature_at_sea_floor';
long_name          = 'Sea Water Potential Temperature at Sea Floor';
units              = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
%% 
ncid = netcdf.open([fpath 'ukesm1-0-ll_r1i1p1f2_historical_tob_60arcmin_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
tob(tob >= 1.00e+20) = NaN;


%% Time
yr = ((time)/12)+1601;

ztest = squeeze(tob(:,:,230));

figure
pcolor(ztest)
shading flat

%%
temp_btm = double(tob);

clear tob

save([fpath 'ukesm_hist_temp_btm_monthly_1850_2014.mat'],'temp_btm','yr',...
    'long_name','standard_name','units','lat','lon','time');





