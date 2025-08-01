% Read CMIP6 netcdfs
% UK-ESM SSP 534-over DOES NOT HAVE ZMESO
% NetCDF on the DKRZ ISIMIP server
% Int top 200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';

%% Meso Zoop zall
ncdisp([fpath 'ukesm1-0-ll_r4i1p1f2_ssp534-over_zmeso_60arcmin_global_monthly_2101_2300.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1' 
%calendar   = '360_day'
                         
%%
ncid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp534-over_zmeso_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Vars
%lat: 330
%lon: 360
%lev: 75
%zmeso: 360x180x75x1032
%time: 1980
%NaNs = 1.000000020040877e+20

%% Get subsets of zmeso so smaller memory
% Time
yr = ((time+1)/12)+1601;
%z200 = find(lev <= 200);
load([fpath 'ukesm_ssp534_zlev.mat'],'z200','zlev');

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,0],[360 180 length(z200) length(time)]);
    
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'ukesm1-0-ll_r4i1p1f2_ssp534-over_thkcello_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,0,0],[360 180 length(z200) length(time)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Make sure oriented the same
whos zmeso thkcello
ztest = squeeze(zmeso(:,:,1,520));
ttest = squeeze(thkcello(:,:,1,520));

figure
pcolor(ztest)
figure
pcolor(ttest)

%% both flipped L-R from what it should be
% zmeso = fliplr(zmeso);
% thkcello = fliplr(thkcello);

%% Integrate top 200 m
zmeso_200 = squeeze(sum((zmeso.*thkcello),3,'omitnan'));

%%
whos zmeso_200
ztest2 = squeeze(zmeso_200(:,:,700));

figure
pcolor(ztest2)
shading flat

%%
%clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ukesm_ssp534-over_zmeso_200_monthly_2101_2300.mat'],'zmeso_200',...
    'yr','time','long_name','standard_name','lat','lon',...
    'units_orig','units_vint','z200','zlev');






