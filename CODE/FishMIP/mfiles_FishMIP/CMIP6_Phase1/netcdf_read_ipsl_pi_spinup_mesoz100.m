% Read Fish-MIP netcdfs
% IPSL PreIndust
% Re-do vertical integration of top 100m
% Need to account for thkcello (cell thickness)
% Before was even 10 m top 100, but not anymore

clear all
close all

spath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/IPSL/preindust/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/preindust/';

%% Meso Zoop zall
ncdisp([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_zmeso_onedeg_global_monthly_1601_2100.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton expressed as Carbon in sea water';
units         = 'mol m-3';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_zmeso_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subset
% Time
%spin 1850-1949, hist 1950-2014, and ssp 2015-2100
%1200 mo, 780 mo, 1032 mo 
yr = ((time)/12)+1601; %months since 1601-1-1 = first month is Jan 1601

spin = find(yr>=1850 & yr<1950);
runs = find(yr>=1950);

runs(end)
runs(1)
spin(end)

% Top 100 m
z100 = find(olevel <= 100);

%%
i = nvars;
zmeso = netcdf.getVar(ncid,i-1, [0,0,0,spin(1)-1],[360 180 length(z100) length(spin)]);
zmeso(zmeso >= 1.00e+20) = NaN;
netcdf.close(ncid);

% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_thkcello_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);
% just last var = thkcello
t = tvars;
thkcello = netcdf.getVar(tcid,t-1, [0,0,0,spin(1)-1],[360 180 length(z100) length(spin)]);
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;

% check orientation
test=squeeze(double(zmeso(:,:,350)));
figure
pcolor(test)
colorbar

test1=squeeze(double(thkcello(:,:,350)));
figure
pcolor(test1)
colorbar

%% Integrate top 200 m
zmeso_100 = squeeze(nansum((zmeso.*thkcello),3));

test2=squeeze(double(zmeso_100(:,:,350)));
pcolor(test2)
colorbar

%%
zmeso_100 = fliplr(zmeso_100);
test=squeeze(double(zmeso_100(:,:,200)));
pcolor(test)
colorbar

%%
clear zmeso thkcello

units_orig = units;
units_vint = 'mol m-2';

%%
save([fpath 'ipsl_pi_zmeso_100_monthly_1850_1949.mat'],'zmeso_100','time',...
    'yr','spin','long_name','standard_name','units_orig','units_vint','olevel','z100');
save([spath 'ipsl_pi_zmeso_100_monthly_1850_1949.mat'],'zmeso_100','time',...
    'yr','spin','long_name','standard_name','units_orig','units_vint','olevel','z100');





