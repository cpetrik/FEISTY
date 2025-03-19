% Read CMIP6 netcdfs
% IPSL SSP 585
% Vert int top 200 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';

% Meso Zoop zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_zmeso_60arcmin_global_monthly_2101_2300.nc'])

%
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton Expressed as Carbon in Sea Water';
units         = 'mol m-3';
missing_value = 1.000000020040877e+20;
%Size:       360x180x75x1032
%Dimensions: i,j,time
%time units = 'months since 1601-1-1'
                         
%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_zmeso_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

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
    zmeso = netcdf.getVar(ncid,n-1, [0,0,0,0], [360 180 length(z200) length(time)]);
end
netcdf.close(ncid);
zmeso(zmeso >= 1.00e+20) = NaN;

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
whos zmeso thkcello
ztest = squeeze(zmeso(:,:,1,150));
ttest = squeeze(thkcello(:,:,1,150));

figure
pcolor(ztest); shading flat
figure
pcolor(ttest); shading flat

%% both flipped L-R from what it should be
% zmeso = fliplr(zmeso);
% thkcello = fliplr(thkcello);

%% Integrate top 200 m
zmeso_200 = squeeze(sum((zmeso.*thkcello),3,'omitnan'));

%%
whos zmeso_200
ztest2 = squeeze(zmeso_200(:,:,100));

figure
pcolor(ztest2)
shading flat

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

save([fpath 'ipsl_ssp585_zmeso_200_monthly_2101_2300.mat'],'zmeso_200','yr',...
    'long_name','standard_name','units_orig','units_vint','lat','lon',...
    'time','z200','olevel');





