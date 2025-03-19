% Plot ssp126

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';

%%
load([fpath 'Means_ipsl_ssp126_monthly_2015_2100.mat']);
load([fpath 'Means_ipsl_ssp126_monthly_2101_2300.mat']);

%%
figure
subplot(2,2,1)
plot(ssp126_yr1,ssp126_Tp1,'r'); hold on
plot(ssp126_yr2,ssp126_Tp2,'r'); hold on

subplot(2,2,2)
plot(ssp126_yr1,ssp126_Tb1,'b'); hold on
plot(ssp126_yr2,ssp126_Tb2,'b'); hold on

subplot(2,2,3)
plot(ssp126_yr1,ssp126_Zm1,'color',[0.75 0 0.5]); hold on
plot(ssp126_yr2,ssp126_Zm2,'color',[0.75 0 0.5]); hold on

subplot(2,2,4)
plot(ssp126_yr1,ssp126_Det1,'color',[0 0.5 0.75]); hold on
plot(ssp126_yr2,ssp126_Det2,'color',[0 0.5 0.75]); hold on

