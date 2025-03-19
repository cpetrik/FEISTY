% Plot all scenarios together

clear 
close all

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% Hist
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

load([hpath 'Means_ukesm_hist_monthly.mat']);

%% SSP 126
bpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp126/';

load([bpath 'Means_ukesm_ssp126_monthly_2015_2150.mat']);
load([bpath 'Means_ukesm_ssp126_monthly_2151_2300.mat']);

ssp126_yr = [ssp126_yr1; ssp126_yr2];

ssp126_Tp = [ssp126_Tp1 ssp126_Tp2];

ssp126_Tb = [ssp126_Tb1 ssp126_Tb2];

ssp126_Zm = [ssp126_Zm1 ssp126_Zm2];

ssp126_Det = [ssp126_Det1 ssp126_Det2];

%% SSP 585
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';

load([wpath 'Means_ukesm_ssp585_monthly_2015_2150.mat']);
load([wpath 'Means_ukesm_ssp585_monthly_2151_2300.mat']);

ssp585_yr = [ssp585_yr1; ssp585_yr2];

ssp585_Tp = [ssp585_Tp1 ssp585_Tp2];

ssp585_Tb = [ssp585_Tb1 ssp585_Tb2];

ssp585_Zm = [ssp585_Zm1 ssp585_Zm2];

ssp585_Det = [ssp585_Det1 ssp585_Det2];

%% SSP 534
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';

load([wpath 'Means_ukesm_ssp534_monthly_2040_2100.mat']);

%% Rolling means ~ annual
hist_yr = movmean(hist_yr,12);
ssp126_yr = movmean(ssp126_yr,12);
ssp585_yr = movmean(ssp585_yr,12);
ssp534_yr = movmean(ssp534_yr,12);

hist_Tp = movmean(hist_Tp,12);
ssp126_Tp = movmean(ssp126_Tp,12);
ssp585_Tp = movmean(ssp585_Tp,12);
ssp534_Tp = movmean(ssp534_Tp,12);

hist_Tb = movmean(hist_Tb,12);
ssp126_Tb = movmean(ssp126_Tb,12);
ssp585_Tb = movmean(ssp585_Tb,12);
ssp534_Tb = movmean(ssp534_Tb,12);

hist_Zm = movmean(hist_Zm,12);
ssp126_Zm = movmean(ssp126_Zm,12);
ssp585_Zm = movmean(ssp585_Zm,12);
%ssp534_Zm = movmean(ssp534_Zm,12);

hist_Det = movmean(hist_Det,12);
ssp126_Det = movmean(ssp126_Det,12);
ssp585_Det = movmean(ssp585_Det,12);
ssp534_Det = movmean(ssp534_Det,12);

%%
figure
subplot(2,2,1)
plot(hist_yr,hist_Tp,'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tp,'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tp,'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tp,'color',[0 0.75 0.5],'LineWidth',2);
title('UKESM T pel')
legend('Hist','126','585','534')
legend('location','northwest')

subplot(2,2,2)
plot(hist_yr,hist_Tb,'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tb,'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tb,'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tb,'color',[0 0.75 0.5],'LineWidth',2);
title('UKESM T btm')

% subplot(2,2,3)
% plot(hist_yr,hist_Zm,'k','LineWidth',1.5); hold on
% plot(ssp126_yr,ssp126_Zm,'b','LineWidth',1.5); hold on
% plot(ssp585_yr,ssp585_Zm,'r','LineWidth',1.5); hold on
% plot(ssp534_yr,ssp534_Zm,'color',[0 0.75 0.5],'LineWidth',1.5);
% title('UKESM Zmeso')

subplot(2,2,4)
plot(hist_yr,hist_Det,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Det,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Det,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Det,'color',[0 0.75 0.5],'LineWidth',1.5);
title('UKESM Det btm')
print('-dpng',[pp 'UKESM_forcing_global_means_all_scenarios.png'])