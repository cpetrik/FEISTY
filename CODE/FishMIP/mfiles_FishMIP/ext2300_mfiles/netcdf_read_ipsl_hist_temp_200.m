% Read CMIP6 netcdfs
% IPSL Hist
% Vert mean top 200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/hist/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/hist/';

%% temp zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thetao_onedeg_global_monthly_1850_2014.nc'])

%%
standard_name      = 'sea_water_potential_temperature';
long_name          = 'Sea Water Potential Temperature';
units              = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
%% 
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thetao_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end


%% Time
yr = ((time)/12)+1601;
%runs = find(yr>1950 & yr<=2015);
z200 = find(olevel <= 200);

%% Get subset of time & depth
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    thetao = netcdf.getVar(ncid,n-1, [0,0,0,0], [360 180 length(z200) 1980]);
end
netcdf.close(ncid);
thetao(thetao >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_historical_thkcello_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[360 180 length(z200) 1980]);
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

%% Depth-weighted mean top 200 m
temp_200 = squeeze( sum((thetao.*thkcello),3,'omitnan') ./ sum((thkcello),3,'omitnan') );

%%
whos temp_200
ztest2 = squeeze(temp_200(:,:,100));

figure
pcolor(ztest2)
shading flat

%%
clear thetao thkcello

units_orig = units;

save([fpath 'ipsl_hist_temp_200_monthly_1850_2014.mat'],'temp_200','yr',...
    'long_name','standard_name','units_orig','units','lat','lon',...
    'z200','olevel');





