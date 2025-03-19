% Plot ssp585

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';

%%
load([fpath 'Means_ipsl_ssp585_monthly_2015_2100.mat']);
load([fpath 'Means_ipsl_ssp585_monthly_2101_2300.mat']);

%%
figure
subplot(2,2,1)
plot(ssp585_yr1,ssp585_Tp1,'r'); hold on
plot(ssp585_yr2,ssp585_Tp2,'r'); hold on

subplot(2,2,2)
plot(ssp585_yr1,ssp585_Tb1,'b'); hold on
plot(ssp585_yr2,ssp585_Tb2,'b'); hold on

subplot(2,2,3)
plot(ssp585_yr1,ssp585_Zm1,'color',[0.75 0 0.5]); hold on
plot(ssp585_yr2,ssp585_Zm2,'color',[0.75 0 0.5]); hold on

subplot(2,2,4)
plot(ssp585_yr1,ssp585_Det1,'color',[0 0.5 0.75]); hold on
plot(ssp585_yr2,ssp585_Det2,'color',[0 0.5 0.75]); hold on

