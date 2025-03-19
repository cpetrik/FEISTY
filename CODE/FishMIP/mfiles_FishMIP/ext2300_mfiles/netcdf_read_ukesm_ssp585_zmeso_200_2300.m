% Read CMIP6 netcdfs
% UK-ESM SSP 585
% Mean top 200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';

%% Temp zall
ncdisp([fpath 'ukesm1-0-ll_r4i1p1f2_ssp585_zmeso_60arcmin_global_monthly_2015_2300.nc'])

%%
standard_name  = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name      = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units          = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x3432
%Dimensions: i,j,time
%time units = 'months since 1601-1-1' 
%calendar   = '360_day'
                         
%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp585_zmeso_60arcmin_global_monthly_2015_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of zmeso so smaller memory 
%360x180x31x1811
% Time
yr = ((time)/12)+1601;
runs = find(yr>2150);
z200 = find(lev <= 200);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
    
end
netcdf.close(ncid);

zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp585_thkcello_60arcmin_global_monthly_2015_2300.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% check if either are L-R flipped (lat,lon) 
whos thkcello zmeso
ttest = squeeze(thkcello(:,:,1,100));
ztest = squeeze(zmeso(:,:,1,100));

figure
pcolor(ztest)
shading flat

figure
pcolor(ttest)
shading flat

%% Integrate top 200 m
zmeso_200 = squeeze(sum((zmeso.*thkcello),3,'omitnan')) ;

%%
whos zmeso_200
ztest2 = squeeze(zmeso_200(:,:,1001));

figure
pcolor(ztest2)
shading flat

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ukesm_ssp585_zmeso_200_monthly_2150_2300.mat'],'zmeso_200',...
    'yr','time','long_name','standard_name','lat','lon',...
    'units_orig','units_vint','z200','lev','runs');

