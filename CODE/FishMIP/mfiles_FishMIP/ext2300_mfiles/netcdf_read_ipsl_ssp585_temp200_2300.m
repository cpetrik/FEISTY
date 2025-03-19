% Read CMIP6 netcdfs
% IPSL SSP 585
% Vert mean top 200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';

% Temp zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_thetao_60arcmin_global_monthly_2101_2300.nc'])

%
	standard_name      = 'sea_water_potential_temperature';
	long_name          = 'Sea Water Potential Temperature';
	units              = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
% 
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_thetao_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


% Time
yr = ((time+1)/12)+1601;
z200 = find(olevel <= 200);

% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    thetao = netcdf.getVar(ncid,n-1, [0,0,0,0], [360 180 length(z200) length(time)]);
end
netcdf.close(ncid);
thetao(thetao >= 1.00e+20) = NaN;

% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_thkcello_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[360 180 length(z200) length(time)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

% Make sure oriented the same
whos thetao thkcello
ztest = squeeze(thetao(:,:,1,150));
ttest = squeeze(thkcello(:,:,1,150));

figure
pcolor(ztest); shading flat
figure
pcolor(ttest); shading flat


%% Integrate top 200 m
temp_200 = squeeze( sum((thetao.*thkcello),3,'omitnan') ./ sum((thkcello),3,'omitnan') );

%
whos temp_200
ztest2 = squeeze(temp_200(:,:,100));

figure
pcolor(ztest2)
shading flat

%%
clear thetao thkcello

units_orig = units;

save([fpath 'ipsl_ssp585_temp_200_monthly_2101_2300.mat'],'temp_200','yr',...
    'long_name','standard_name','units_orig','units','lat','lon',...
    'time','z200','olevel');





