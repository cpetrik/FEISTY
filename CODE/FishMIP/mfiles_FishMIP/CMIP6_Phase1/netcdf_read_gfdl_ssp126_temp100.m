% Read Fish-MIP netcdfs
% GFDL ssp126
% Re-do mean of top 100m to be depth-weighted
% Need to account for thkcello (cell thickness)
% Before was even 10 m top 100, but not anymore

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/ssp126/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp126/';

%% thetao
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_ssp126_thetao_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp126_thetao_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subset of temp
% Time
yr = ((time)/12)+1601-1;
z100 = find(lev <= 100);

i = nvars;
thetao = netcdf.getVar(ncid,i-1, [0,0,0,0],[360 180 length(z100) length(time)]);
thetao(thetao >= 1.00e+20) = NaN;
netcdf.close(ncid);

%% Subset of thkcello
tcid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_ssp126_thkcello_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
t = tvars;
thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[360 180 length(z100) length(time)]);
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Mean top 100 m
temp_100 = squeeze(nansum((thetao.*thkcello),3)) ./ squeeze(nansum(thkcello,3));

test=squeeze(double(temp_100(:,:,60)));
pcolor(test)
colorbar

%% Time

save([fpath 'gfdl_ssp126_temp_100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','lev','z100');
save([spath 'gfdl_ssp126_temp_100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','lev','z100');





