% Plot all scenarios together

clear 
close all

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% Hist
hpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';

load([hpath 'Means_cesm2_hist_zmeso_monthly_1850_2014.mat']);

%% SSP 126
bpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp126/';

load([bpath 'Means_cesm2_ssp126_zmeso_monthly_2015_2299.mat']);

%% SSP 585
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp585/';

load([wpath 'Means_cesm2_ssp585_zmeso_monthly_2015_2299.mat']);

%% SSP 534
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';

load([wpath 'Means_cesm2_ssp534_zmeso_monthly_2040_2299.mat']);

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
title('CESM T pel')
legend('Hist','126','585','534')
legend('location','northwest')

subplot(2,2,2)
plot(hist_yr,hist_Tb,'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tb,'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tb,'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tb,'color',[0 0.75 0.5],'LineWidth',2);
title('CESM T btm')

subplot(2,2,3)
plot(hist_yr,hist_Zm,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Zm,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Zm,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Zm,'color',[0 0.75 0.5],'LineWidth',1.5);
title('CESM Zmeso')

subplot(2,2,4)
plot(hist_yr,hist_Det,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Det,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Det,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Det,'color',[0 0.75 0.5],'LineWidth',1.5);
title('CESM Det btm')
print('-dpng',[pp 'CESM2-WACCM_forcing_zmeso_global_means_all_scenarios.png'])