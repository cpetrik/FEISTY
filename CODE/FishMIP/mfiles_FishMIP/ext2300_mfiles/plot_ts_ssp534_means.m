% Plot 534

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';

%%
load([fpath 'Means_ipsl_ssp534_monthly_2040_2100.mat']);

Tp1 = ssp534_Tp;
Tb1 = ssp534_Tb;
Zm1 = ssp534_Zm;
Det1 = ssp534_Det;

clear  ssp534_Tp ssp534_Tb ssp534_Zm ssp534_Det 

%%
load([fpath 'Means_ipsl_ssp534_monthly_2101_2300.mat']);

%%
figure
subplot(2,2,1)
plot(ssp534_yr,Tp1,'r'); hold on
plot(ssp534_yr2,ssp534_Tp,'r'); hold on

subplot(2,2,2)
plot(ssp534_yr,Tb1,'b'); hold on
plot(ssp534_yr2,ssp534_Tb,'b'); hold on

subplot(2,2,3)
plot(ssp534_yr,Zm1,'color',[0.75 0 0.5]); hold on
plot(ssp534_yr2,ssp534_Zm,'color',[0.75 0 0.5]); hold on

subplot(2,2,4)
plot(ssp534_yr,Det1,'color',[0 0.5 0.75]); hold on
plot(ssp534_yr2,ssp534_Det,'color',[0 0.5 0.75]); hold on

%%
load([fpath 'ipsl_ssp534-over_det_monthly_2040_2300.mat'],'expc')

det = double(expc);
test4 = squeeze(det(:,:,70));
% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));

[ni,nj,nt] = size(det);

Det= double(reshape(det,ni*nj,nt));
Det= Det(WID,:);
Det2 = mean(Det);

%%
figure
plot(Det2,'color',[0 0.5 0.75]); hold on