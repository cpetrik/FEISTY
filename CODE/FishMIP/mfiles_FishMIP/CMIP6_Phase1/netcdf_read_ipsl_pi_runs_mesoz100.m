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

%zoo: 360x180x75x6000

%% Time
%spin 1850-1949, hist 1950-2014, and ssp 2015-2100
%1200 mo, 780 mo, 1032 mo 
yr = ((time)/12)+1601; %months since 1601-1-1 = first month is Jan 1601

spin = find(yr>=1850 & yr<1950);
runs = find(yr>=1950);

runs(end)
runs(1)
spin(end)

z100 = find(olevel <= 100);

%split into smaller chunks
ct = length(runs)/3;
run1 = runs(1:ct);
run2 = runs(ct+1:2*ct);
run3 = runs(2*ct+1:end);

%% Get subset of mesoz & thkcello 
i = nvars;
zmeso1 = netcdf.getVar(ncid,i-1, [0,0,0,run1(1)-1],[360 180 length(z100) length(run1)]);
zmeso1(zmeso1 >= 1.00e+20) = NaN;

% Subset of thkcello
tcid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_picontrol_thkcello_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);
%just last var = thkcello
t = tvars;
thkcello1 = netcdf.getVar(tcid,t-1, [0,0,0,run1(1)-1],[360 180 length(z100) length(run1)]);
thkcello1(thkcello1 >= 1.00e+20) = NaN;

% check orientation
test=squeeze(double(zmeso1(:,:,250)));
figure
pcolor(test)
colorbar

test1=squeeze(double(thkcello1(:,:,250)));
figure
pcolor(test1)
colorbar

% Mean top 100 m 
zmeso_1 = squeeze(nansum((zmeso1.*thkcello1),3));

%%
clear zmeso1 thkcello1

%% Get subset of zmeso & thkcello
i = nvars;
zmeso2 = netcdf.getVar(ncid,i-1, [0,0,0,run2(1)-1],[360 180 length(z100) length(run2)]);
zmeso2(zmeso2 >= 1.00e+20) = NaN;

% Subset of thkcello
%just last var = thkcello
t = tvars;
thkcello2 = netcdf.getVar(tcid,t-1, [0,0,0,run2(1)-1],[360 180 length(z100) length(run2)]);
thkcello2(thkcello2 >= 1.00e+20) = NaN;

% check orientation
test=squeeze(double(zmeso2(:,:,50)));
figure
pcolor(test)
colorbar

test1=squeeze(double(thkcello2(:,:,50)));
figure
pcolor(test1)
colorbar

% Mean top 100 m 
zmeso_2 = squeeze(nansum((zmeso2.*thkcello2),3));

%%
clear zmeso2 thkcello2

%% Get subset of zmeso & thkcello
i = nvars;
zmeso3 = netcdf.getVar(ncid,i-1, [0,0,0,run3(1)-1],[360 180 length(z100) length(run3)]);
zmeso3(zmeso3 >= 1.00e+20) = NaN;
netcdf.close(ncid);

% Subset of thkcello
%just last var = thkcello
t = tvars;
thkcello3 = netcdf.getVar(tcid,t-1, [0,0,0,run3(1)-1],[360 180 length(z100) length(run3)]);
thkcello3(thkcello3 >= 1.00e+20) = NaN;
netcdf.close(tcid);

% check orientation
test=squeeze(double(zmeso3(:,:,20)));
figure
pcolor(test)
colorbar

test1=squeeze(double(thkcello3(:,:,20)));
figure
pcolor(test1)
colorbar

% Mean top 100 m 
zmeso_3 = squeeze(nansum((zmeso3.*thkcello3),3));

%%
clear zmeso3 thkcello3

%% put together
zmeso_10 = cat(3, zmeso_1, zmeso_2);
zmeso_100 = cat(3, zmeso_10, zmeso_3);

%%
zmeso_100 = fliplr(zmeso_100);
test2=squeeze(double(zmeso_100(:,:,250)));
figure
pcolor(test2)
colorbar

%%
units_orig = units;
units_vint = 'mol m-2';

%%
save([fpath 'ipsl_pi_zmeso_100_monthly_1950_2100.mat'],'zmeso_100','time',...
    'yr','runs','long_name','standard_name','units_orig','units_vint','olevel','z100');
save([spath 'ipsl_pi_zmeso_100_monthly_1950_2100.mat'],'zmeso_100','time',...
    'yr','runs','long_name','standard_name','units_orig','units_vint','olevel','z100');





