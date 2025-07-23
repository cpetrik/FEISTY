% Make mat files of interpolated time series from CESM
% Hist 1950-2014
% 200m vertical integrations

clear 
close all

wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% bias corrected mesozoo from zooc and diatfrac

load([fpath 'cesm2_hist_biascorr_zooc_zmeso_diatfraconly_molC_monthly_1850_2014.mat'])

load([wpath 'Data_grid_cesm2_cmip6_2300.mat']);

%%  Means over all grid cells
[ni,nj,nt] = size(zmeso_corr);
Zm1 = double(reshape(zooc_corr,ni*nj,nt));
Zm2 = double(reshape(zmeso_corr,ni*nj,nt));

Zm1 = Zm1(GRD.ID,:);
Zm2 = Zm2(GRD.ID,:);

hist_Zm1 = mean(Zm1,'omitnan');
hist_Zm2 = mean(Zm2,'omitnan');

%% Plot t.s. with other ESMs

Chist_yr = movmean(yr,12);
Chist_Zm1 = movmean(hist_Zm1,12);
Chist_Zm2 = movmean(hist_Zm2,12);

%% Hist UKESM
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

load([hpath 'Means_ukesm_hist_monthly.mat'],'hist_yr','hist_Zm');

Uhist_yr = movmean(hist_yr,12);
Uhist_Zm = movmean(hist_Zm,12);

clear hist_yr hist_Zm

%% Hist IPSL
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/hist/';

load([hpath 'Means_ipsl_hist_monthly.mat'],'hist_yr','hist_Zm');

Ihist_yr = movmean(hist_yr,12);
Ihist_Zm = movmean(hist_Zm,12);

clear hist_yr hist_Zm

%%
figure
plot(Uhist_yr,Uhist_Zm,'b','LineWidth',2); hold on
plot(Ihist_yr,Ihist_Zm,'r','LineWidth',2); hold on
plot(Chist_yr,Chist_Zm1,'color',[0 0.75 0.5],'LineWidth',2);
plot(Chist_yr,Chist_Zm2,'color',[0.75 0 0.5],'LineWidth',2);
title('Mesozoo')
legend('UKESM','IPSL','bcCESMzooc','bcCESMzmeso')
legend('location','northwest')
print('-dpng',[ppath 'ts_ESMs_hist_biascorr_zoo_test.png'])














