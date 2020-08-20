% Read Fish-MIP netcdfs
% GFDL PreIndust

clear all
close all

fpath='/Users/cpetrik/Dropbox/ESM_data/Fish-MIP/CMIP6/GFDL/preindust/';

%% Meso Zoop zall
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_picontrol_zmeso_onedeg_global_monthly_1601_2100.nc'])
standard_name = 'mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water';
long_name     = 'Mole Concentration of Mesozooplankton expressed as Carbon in sea water';
units         = 'mol m-3';

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_picontrol_zmeso_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Vars
%lat: 180
%lon: 360
%lev: 35
%zoo: 360x180x35x6000
%time: 6000
%NaNs = 1.0000e+20

%% Get subset of chl
% Time
yr = ((time+1)/12)+1601-1;
spin = find(yr>1850 & yr<=1950);
btwn = find(yr>1900 & yr<=2000);
runs = find(yr>1950 & yr<=2100);
z100 = find(lev <= 100);

i = nvars;
zmeso1 = netcdf.getVar(ncid,i-1, [0,0,0,spin(1)-1],[360 180 length(z100) length(spin)]);
zmeso2 = netcdf.getVar(ncid,i-1, [0,0,0,btwn(1)-1],[360 180 length(z100) length(btwn)]);
zmeso3 = netcdf.getVar(ncid,i-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
%zmeso(zmeso >= 1.00e+20) = NaN;
netcdf.close(ncid);

%%
zm1 = squeeze(zmeso1(:,:,6,:));
zm2 = squeeze(zmeso2(:,:,6,:));
zm3 = squeeze(zmeso3(:,:,6,:));

%% Mean top 100 m 
zm1(zm1 >= 1.00e+20) = NaN;
zm2(zm2 >= 1.00e+20) = NaN;
zm3(zm3 >= 1.00e+20) = NaN;

% zmeso1_100 = squeeze(nanmean(zmeso1,3));
% zmeso2_100 = squeeze(nanmean(zmeso2,3));
% zmeso3_100 = squeeze(nanmean(zmeso3,3));

zmeso1_100 = squeeze(nanmean(zm1,3));
zmeso2_100 = squeeze(nanmean(zm2,3));
zmeso3_100 = squeeze(nanmean(zm3,3));

%%
[ni,nj,nt1] = size(zm1);
[ni,nj,nt2] = size(zm2);
[ni,nj,nt3] = size(zm3);
mz1 = reshape(zm1,ni*nj,nt1);
mz2 = reshape(zm2,ni*nj,nt2);
mz3 = reshape(zm3,ni*nj,nt3);
mmz1 = nanmean(mz1);
mmz2 = nanmean(mz2);
mmz3 = nanmean(mz3);

figure
plot(yr(spin),mmz1,'k'); hold on;
plot(yr(btwn),mmz2,'b'); hold on;
plot(yr(runs),mmz3,'r'); hold on;

%%
figure
pcolor(zmeso1_100')
shading flat
colorbar
caxis([0 5e-3])
title('spin')

figure
pcolor(zmeso2_100')
shading flat
colorbar
caxis([0 5e-3])
title('btwn')

figure
pcolor(zmeso3_100')
shading flat
colorbar
caxis([0 5e-3])
title('runs')

eq=(zmeso1_100==zmeso2_100);
eq2=(zmeso1_100==zmeso3_100);

%%
zmeso1(zmeso1 >= 1.00e+20) = NaN;
zmeso2(zmeso2 >= 1.00e+20) = NaN;
zmeso3(zmeso3 >= 1.00e+20) = NaN;

zmeso1_100 = squeeze(nanmean(zmeso1,3));
zmeso2_100 = squeeze(nanmean(zmeso2,3));
zmeso3_100 = squeeze(nanmean(zmeso3,3));

%%
nz1 = reshape(zmeso1_100,ni*nj,nt1);
nz2 = reshape(zmeso2_100,ni*nj,nt2);
nz3 = reshape(zmeso3_100,ni*nj,nt3);
mnz1 = nanmean(nz1);
mnz2 = nanmean(nz2);
mnz3 = nanmean(nz3);

figure
plot(yr(spin),mnz1,'k'); hold on;
plot(yr(btwn),mnz2,'b'); hold on;
plot(yr(runs),mnz3,'r'); hold on;

%%
save([fpath 'gfdl_pi_zmeso100_monthly_1850_1950_2100.mat'],'time',...
    'yr','runs','long_name','standard_name','units','lev','z100',...
    'zmeso1_100','zmeso2_100','zmeso3_100','spin','btwn');


