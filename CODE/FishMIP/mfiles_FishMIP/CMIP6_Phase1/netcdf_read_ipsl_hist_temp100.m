% Read Fish-MIP netcdfs
% IPSL Hist temp
% Re-do mean of top 100m to be depth-weighted
% Need to account for thkcello (cell thickness)
% Before was even 10 m top 100, but not anymore

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';

%% temp
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thetao_onedeg_global_monthly_1850_2014.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thetao_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subset of thetao
% Time
%yr = ((time+1)/12)+1601-1;
yr = ((time)/12)+1601-1;
runs = find(yr>1949 & yr<=2015);
z100 = find(olevel <= 100);

i = nvars;
thetao = netcdf.getVar(ncid,i-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
thetao(thetao >= 1.00e+20) = NaN;
netcdf.close(ncid);

%% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thkcello_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
t = tvars;
thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Mean top 100 m 
temp_100 = squeeze(nansum((thetao.*thkcello),3)) ./ squeeze(nansum(thkcello,3));

test=squeeze(double(temp_100(:,:,300)));
pcolor(test)
colorbar

%%
save([fpath 'ipsl_hist_temp_100_monthly_1950_2014.mat'],'temp_100','time',...
    'yr','runs','long_name','standard_name','units','olevel','z100');
save([spath 'ipsl_hist_temp_100_monthly_1950_2014.mat'],'temp_100','time',...
    'yr','runs','long_name','standard_name','units','olevel','z100');





