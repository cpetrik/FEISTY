% Read CMIP6 netcdfs
% UK-ESM Hist
% NetCDF on the ESGF server
% In 100-200 m

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

%% zall
ncdisp([fpath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'])

%%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1980
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
%calendar   = '360_day'

%% Time 1
ncid = netcdf.open([fpath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

[ni,nj] = size(longitude);

%% Get subsets of zmeso so smaller memory
% Time
yr1 = ((time+1)/360)+1850;
runs1 = (length(time)-(15*12)+1):length(time);
z100 = find(lev > 100 & lev <= 200);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,z100(1)-1,runs1(1)-1],[ni nj length(z100) length(runs1)]);

end
netcdf.close(ncid);

zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'thkcello_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,z100(1)-1,runs1(1)-1],[ni nj length(z100) length(runs1)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Integrate top 100 m
zmeso_1 = squeeze(nansum((zmeso.*thkcello),3));

%%
clear zmeso thkcello


%% Time 2
ncid = netcdf.open([fpath 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of zmeso so smaller memory
% Time
yr2 = ((time+1)/360)+1850;
runs2 = 1:length(time);
z100 = find(lev > 100 & lev <= 200);

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1, [0,0,z100(1)-1,runs2(1)-1],[ni nj length(z100) length(runs2)]);

end
netcdf.close(ncid);

zmeso(zmeso >= 1.00e+20) = NaN;

%% Subset of thkcello
tcid = netcdf.open([fpath 'thkcello_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_200001-201412.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    thkcello = netcdf.getVar(tcid,t-1, [0,0,z100(1)-1,runs2(1)-1],[ni nj length(z100) length(runs2)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

%% Integrate
zmeso_2 = squeeze(nansum((zmeso.*thkcello),3));

zmeso_200 = cat(3,zmeso_1,zmeso_2);
yr = [yr1;yr2];
runs = [runs1,runs1(end)+runs2];

%%
clear zmeso thkcello zmeso_1 zmeso_2

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ukesm_hist_zmeso_100-200_monthly_1985_2014.mat'],'zmeso_200',...
    'yr1','yr2','yr','runs','long_name','standard_name','latitude','longitude',...
    'units_orig','units_vint','z100','lev');
