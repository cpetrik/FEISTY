% Make mat files of interpolated time series from CESM
% Hist 1950-2014
% 200m vertical integrations

clear 
close all

gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% bias corrected mesozoo from zooc and diatfrac
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';

load([fpath 'cesm2_hist_biascorr_zooc_zmeso_diatfraconly_molC_monthly_1850_2014.mat'])
load([fpath 'Means_cesm2_hist_monthly_1850_2014.mat'],'hist_yr');

load([gpath 'Data_grid_cesm2_cmip6_2300.mat']);

hist_zooc = zooc_corr;
hist_zmeso = zmeso_corr;

%% SSP 126
bpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp126/';

load([bpath 'cesm2_ssp126_zooc_150_monthly_2015_2299.mat']);
load([bpath 'cesm2_ssp126_zmeso_150_monthly_2015_2299.mat']);
load([bpath 'Means_cesm2_ssp126_monthly_2015_2299.mat'],'ssp126_yr2');

[ni,nj,nt] = size(zmeso_150);

s126_zooc = (zooc_150) - repmat(diffZooc,1,1,nt);
s126_zmeso = (zmeso_150) - repmat(diffZmeso,1,1,nt);

%% SSP 585
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp585/';

load([wpath 'cesm2_ssp585_zooc_150_monthly_2015_2299.mat']);
load([wpath 'cesm2_ssp585_zmeso_150_monthly_2015_2299.mat']);
load([wpath 'Means_cesm_ssp585_monthly_2015_2299.mat'],'ssp585_yr1');

[ni,nj,nt] = size(zmeso_150);

s585_zooc = (zooc_150) - repmat(diffZooc,1,1,nt);
s585_zmeso = (zmeso_150) - repmat(diffZmeso,1,1,nt);

%% SSP 534
spath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';

load([spath 'cesm2_ssp534-over_zooc_150_monthly_2040_2299.mat']);
load([spath 'cesm2_ssp534-over_zmeso_150_monthly_2040_2299.mat']);
load([spath 'Means_cesm2_ssp534_monthly_2040_2299.mat'],'ssp534_yr');

[ni,nj,nt] = size(zmeso_150);

s534_zooc = (zooc_150) - repmat(diffZooc,1,1,nt);
s534_zmeso = (zmeso_150) - repmat(diffZmeso,1,1,nt);

%%  Means over all grid cells
Zz1 = double(reshape(hist_zooc,ni*nj,1980));
Zm1 = double(reshape(hist_zmeso,ni*nj,1980));

Zz2 = double(reshape(s126_zooc,ni*nj,3420));
Zm2 = double(reshape(s126_zmeso,ni*nj,3420));

Zz3 = double(reshape(s585_zooc,ni*nj,3420));
Zm3 = double(reshape(s585_zmeso,ni*nj,3420));

Zz4 = double(reshape(s534_zooc,ni*nj,3120));
Zm4 = double(reshape(s534_zmeso,ni*nj,3120));

hist_z = mean(Zz1,'omitnan');
hist_m = mean(Zm1,'omitnan');

s126_z = mean(Zz2,'omitnan');
s126_m = mean(Zm2,'omitnan');

s585_z = mean(Zz3,'omitnan');
s585_m = mean(Zm3,'omitnan');

s534_z = mean(Zz4,'omitnan');
s534_m = mean(Zm4,'omitnan');

%% Rolling means ~ annual
hist_yr = movmean(hist_yr,12);
ssp126_yr = movmean(ssp126_yr2,12);
ssp585_yr = movmean(ssp585_yr1,12);
ssp534_yr = movmean(ssp534_yr,12);

hist_Zz = movmean(hist_z,12);
ssp126_Zz = movmean(s126_z,12);
ssp585_Zz = movmean(s585_z,12);
ssp534_Zz = movmean(s534_z,12);

hist_Zm = movmean(hist_m,12);
ssp126_Zm = movmean(s126_m,12);
ssp585_Zm = movmean(s585_m,12);
ssp534_Zm = movmean(s534_m,12);

%%
figure
subplot(2,2,1)
plot(hist_yr,hist_Zz,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Zz,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Zz,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Zz,'color',[0 0.75 0.5],'LineWidth',1.5);
title('CESM Zooc bias corr')

subplot(2,2,2)
plot(hist_yr,hist_Zm,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Zm,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Zm,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Zm,'color',[0 0.75 0.5],'LineWidth',1.5);
title('CESM Zmeso bias corr')

print('-dpng',[ppath 'CESM2-WACCM_global_zoo_biascorr_all_scenarios.png'])











