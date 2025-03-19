% Plot all scenarios together

clear 
close all

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% Hist
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/hist/';

load([hpath 'Means_ipsl_hist_monthly.mat']);

%% SSP 126
bpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';

load([bpath 'Means_ipsl_ssp126_monthly_2015_2100.mat']);
load([bpath 'Means_ipsl_ssp126_monthly_2101_2300.mat']);

ssp126_yr = [ssp126_yr1; ssp126_yr2];

ssp126_Tp = [ssp126_Tp1 ssp126_Tp2];

ssp126_Tb = [ssp126_Tb1 ssp126_Tb2];

ssp126_Zm = [ssp126_Zm1 ssp126_Zm2];

ssp126_Det = [ssp126_Det1 ssp126_Det2];

%% SSP 585
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';

load([wpath 'Means_ipsl_ssp585_monthly_2015_2100.mat']);
load([wpath 'Means_ipsl_ssp585_monthly_2101_2300.mat']);

ssp585_yr = [ssp585_yr1; ssp585_yr2];

ssp585_Tp = [ssp585_Tp1 ssp585_Tp2];

ssp585_Tb = [ssp585_Tb1 ssp585_Tb2];

ssp585_Zm = [ssp585_Zm1 ssp585_Zm2];

ssp585_Det = [ssp585_Det1 ssp585_Det2];

%% SSP 534
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';

load([wpath 'Means_ipsl_ssp534_monthly_2040_2100.mat']);

ssp534_Tp1 = ssp534_Tp;
ssp534_Tb1 = ssp534_Tb;
ssp534_Zm1 = ssp534_Zm;
ssp534_Det1 = ssp534_Det;

clear  ssp534_Tp ssp534_Tb ssp534_Zm ssp534_Det 

load([wpath 'Means_ipsl_ssp534_monthly_2101_2300.mat']);

ssp534_Tp2 = ssp534_Tp;
ssp534_Tb2 = ssp534_Tb;
ssp534_Zm2 = ssp534_Zm;
ssp534_Det2 = ssp534_Det;

clear  ssp534_Tp ssp534_Tb ssp534_Zm ssp534_Det

ssp534_yr = [ssp534_yr; ssp534_yr2];

ssp534_Tp = [ssp534_Tp1 ssp534_Tp2];

ssp534_Tb = [ssp534_Tb1 ssp534_Tb2];

ssp534_Zm = [ssp534_Zm1 ssp534_Zm2];

ssp534_Det = [ssp534_Det1 ssp534_Det2];

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
ssp534_Zm = movmean(ssp534_Zm,12);

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
title('T pel')
legend('Hist','126','585','534')
legend('location','northwest')

subplot(2,2,2)
plot(hist_yr,hist_Tb,'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tb,'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tb,'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tb,'color',[0 0.75 0.5],'LineWidth',2);
title('T btm')

subplot(2,2,3)
plot(hist_yr,hist_Zm,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Zm,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Zm,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Zm,'color',[0 0.75 0.5],'LineWidth',1.5);
title('Zmeso')

subplot(2,2,4)
plot(hist_yr,hist_Det,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Det,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Det,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Det,'color',[0 0.75 0.5],'LineWidth',1.5);
title('Det btm')
print('-dpng',[pp 'IPSL_forcing_global_means_all_scenarios.png'])