% Read CMIP6 netcdfs
% UK-ESM SSP 126
% Mean top 200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp126/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp126/';

%% Temp zall
ncdisp([fpath 'ukesm1-0-ll_r4i1p1f2_ssp126_thetao_60arcmin_global_monthly_2015_2300.nc'])

%%
standard_name      = 'sea_water_potential_temperature';
long_name          = 'Sea Water Potential Temperature';
units              = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x732
%Dimensions: i,j,time
%time units = 'months since 1601-1-1' 
%calendar   = '360_day'
                         
%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp126_thetao_60arcmin_global_monthly_2015_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of thetao so smaller memory
% Time
yr = ((time)/12)+1601;
%runs1 = find(yr1>1950 & yr1<=2015);
z200 = find(lev <= 200);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    thetao = netcdf.getVar(ncid,n-1, [0,0,0,0],[360 180 length(z200) length(time)]);
    
end
netcdf.close(ncid);

thetao(thetao >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp126_thkcello_60arcmin_global_monthly_2015_2300.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[360 180 length(z200) length(time)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% check if either are L-R flipped (lat,lon) 
whos thkcello thetao
ttest = squeeze(thkcello(:,:,1,100));
ztest = squeeze(thetao(:,:,1,100));

figure
pcolor(ztest)
shading flat

figure
pcolor(ttest)
shading flat

%% Integrate top 200 m
temp_200 = squeeze(sum((thetao.*thkcello),3,'omitnan')) ./ squeeze(sum((thkcello),3,'omitnan'));

%%
whos temp_200
ztest2 = squeeze(temp_200(:,:,100));

figure
pcolor(ztest2)
shading flat

%%
clear thetao thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ukesm_ssp126_temp_200_monthly_2015_2300.mat'],'temp_200',...
    'yr','time','long_name','standard_name','lat','lon',...
    'units','z200','lev');

