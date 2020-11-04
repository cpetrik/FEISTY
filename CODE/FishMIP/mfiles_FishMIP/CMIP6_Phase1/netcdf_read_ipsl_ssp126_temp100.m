% Read Fish-MIP netcdfs
% IPSL ssp126
% Re-do mean of top 100m to be depth-weighted
% Need to account for thkcello (cell thickness)
% Before was even 10 m top 100, but not anymore

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/ssp126/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp126/';

%% thetao
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_thetao_onedeg_global_monthly_2015_2100.nc'])
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_thetao_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subset of thetao
% Time
yr = ((time)/12)+1601-1;
runs = find(yr>2014 & yr<=2100);
z100 = find(olevel <= 100);

i = nvars;
thetao = netcdf.getVar(ncid,i-1, [0,0,0,0],[360 180 length(z100) length(time)]);
thetao(thetao >= 1.00e+20) = NaN;
netcdf.close(ncid);

%% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_thkcello_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
t = tvars;
thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[360 180 length(z100) length(time)]);
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% thkcello is L-R flipped (lat,lon) compared to zmeso
thkcello = fliplr(thkcello);

%% Mean top 100 m 
temp_100 = squeeze(nansum((thetao.*thkcello),3)) ./ squeeze(nansum(thkcello,3));

test=squeeze(double(temp_100(:,:,45)));
pcolor(test)
colorbar

%% save
save([fpath 'ipsl_ssp126_temp_100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','olevel','z100');
save([spath 'ipsl_ssp126_temp_100_monthly_2015_2100.mat'],'temp_100','time',...
    'yr','long_name','standard_name','units','olevel','z100');





