% Read CMIP6 netcdfs
% CESM Hist
% Integrate top 0-200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
tpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';

%% zoo zall
ncdisp([fpath 'cesm2-waccm_r1i1p1f1_historical_zooc_60arcmin_global_monthly_1850_2014.nc'])

%%
standard_name = 'mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Zooplankton Carbon Concentration';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       320x384x15x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1 00:00:00'

%% 
ncid = netcdf.open([fpath 'cesm2-waccm_r1i1p1f1_historical_zooc_60arcmin_global_monthly_1850_2014.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Get subsets of zooc so smaller memory?
yr = ((time+1)/12)+1601;
z200 = find(lev <= 200e2); %lev in cm
%only top 15 levels saved anyway, so works out bring in whole file

%%
zooc(zooc >= 1.00e+20) = NaN;

%% Subset of thkcello #lev_bnds in m (not cm like lev)
[ni,nj,~,nt] = size(zooc);
thkcello = (lev_bnds(2,z200) - lev_bnds(1,z200))';
thkcello = repmat(thkcello,1,ni,nj,nt);
thkcello = permute(thkcello,[2 3 1 4]);

%% Integrate top 150 m
zooc_150 = squeeze(sum((zooc.*thkcello),3,'omitnan'));

%% Clear 
clear zooc thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'cesm2_hist_zooc_150_monthly_1850_2014.mat'],'zooc_150',...
    'yr','long_name','standard_name','units','lat','lon',...
    'lev','lev_bnds','units_orig','units_vint');
