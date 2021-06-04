% Visualize time series CMIP6 inputs for FEISTY
% raw and relative to 1990-2000

clear all
close all

%% Data
ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';
if (~isfolder(ppath))
    mkdir(ppath)
end

ihpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
i2path='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp126/';
i5path='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';

ghpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
g2path='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp126/';
g5path='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';

%% gfdl
load([ghpath 'gfdl_hist_inputs_ts_nanmeans.mat']);
GHist_Tp = mean(Tp_1yr_hist);
GHist_Tb = mean(Tb_1yr_hist);
GHist_det = mean(det_1yr_hist);
GHist_mz = mean(mz_1yr_hist);

load([g2path 'gfdl_ssp126_inputs_ts_nanmeans.mat']);
GSSP126_Tp = mean(Tp_1yr_ssp126);
GSSP126_Tb = mean(Tb_1yr_ssp126);
GSSP126_det = mean(det_1yr_ssp126);
GSSP126_mz = mean(mz_1yr_ssp126);

load([g5path 'gfdl_ssp585_inputs_ts_nanmeans.mat']);
GSSP585_Tp = mean(Tp_1yr_ssp585);
GSSP585_Tb = mean(Tb_1yr_ssp585);
GSSP585_det = mean(det_1yr_ssp585);
GSSP585_mz = mean(mz_1yr_ssp585);

clear Tp_1yr_hist Tb_1yr_hist det_1yr_hist mz_1yr_hist
clear Tp_1yr_ssp126 Tb_1yr_ssp126 det_1yr_ssp126 mz_1yr_ssp126
clear Tp_1yr_ssp585 Tb_1yr_ssp585 det_1yr_ssp585 mz_1yr_ssp585

%% ipsl
load([ihpath 'ipsl_hist_inputs_ts_nanmeans.mat']);
IHist_Tp = mean(Tp_1yr_hist);
IHist_Tb = mean(Tb_1yr_hist);
IHist_det = mean(det_1yr_hist);
IHist_mz = mean(mz_1yr_hist);

load([i2path 'ipsl_ssp126_inputs_ts_nanmeans.mat']);
ISSP126_Tp = mean(Tp_1yr_ssp126);
ISSP126_Tb = mean(Tb_1yr_ssp126);
ISSP126_det = mean(det_1yr_ssp126);
ISSP126_mz = mean(mz_1yr_ssp126);

load([i5path 'ipsl_ssp585_inputs_ts_nanmeans.mat']);
ISSP585_Tp = mean(Tp_1yr_ssp585);
ISSP585_Tb = mean(Tb_1yr_ssp585);
ISSP585_det = mean(det_1yr_ssp585);
ISSP585_mz = mean(mz_1yr_ssp585);

clear Tp_1yr_hist Tb_1yr_hist det_1yr_hist mz_1yr_hist
clear Tp_1yr_ssp126 Tb_1yr_ssp126 det_1yr_ssp126 mz_1yr_ssp126
clear Tp_1yr_ssp585 Tb_1yr_ssp585 det_1yr_ssp585 mz_1yr_ssp585

%% time
yH = hyrs;
yS = s5yrs;

tid = find(yH>1990 & yH<=2000);

GHistTp90 = mean(GHist_Tp(tid));
IHistTp90 = mean(IHist_Tp(tid));

GHistTb90 = mean(GHist_Tb(tid));
IHistTb90 = mean(IHist_Tb(tid));

GHistDet90 = mean(GHist_det(tid));
IHistDet90 = mean(IHist_det(tid));

GHistMz90 = mean(GHist_mz(tid));
IHistMz90 = mean(IHist_mz(tid));

%% Raw change
figure(1)
subplot(2,2,1)
plot(yH,GHist_Tp,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GSSP126_Tp,'b','LineWidth',1); hold on;
plot(yS,GSSP585_Tp,'r','LineWidth',1); hold on;
plot(yH,IHist_Tp,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,ISSP126_Tp,'c','LineWidth',1); hold on;
plot(yS,ISSP585_Tp,'m','LineWidth',1); hold on;
title('Tp')
legend('Ghist','G126','G585','Ihist','I126','I585')
legend('location','northwest')
xlim([1970 2100])

subplot(2,2,2)
plot(yH,GHist_Tb,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GSSP126_Tb,'b','LineWidth',1); hold on;
plot(yS,GSSP585_Tb,'r','LineWidth',1); hold on;
plot(yH,IHist_Tb,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,ISSP126_Tb,'c','LineWidth',1); hold on;
plot(yS,ISSP585_Tb,'m','LineWidth',1); hold on;
title('Tb')
xlim([1970 2100])

subplot(2,2,3)
plot(yH,log10(GHist_mz),'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,log10(GSSP126_mz),'b','LineWidth',1); hold on;
plot(yS,log10(GSSP585_mz),'r','LineWidth',1); hold on;
plot(yH,log10(IHist_mz),'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,log10(ISSP126_mz),'c','LineWidth',1); hold on;
plot(yS,log10(ISSP585_mz),'m','LineWidth',1); hold on;
title('Mesoz')
ylabel('Biomass (log_1_0 g m^-^2)')
xlim([1970 2100])

subplot(2,2,4)
plot(yH,log10(GHist_det),'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,log10(GSSP126_det),'b','LineWidth',1); hold on;
plot(yS,log10(GSSP585_det),'r','LineWidth',1); hold on;
plot(yH,log10(IHist_det),'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,log10(ISSP126_det),'c','LineWidth',1); hold on;
plot(yS,log10(ISSP585_det),'m','LineWidth',1); hold on;
title('Det')
ylabel('Flux (log_1_0 g m^-^2 s^-^1)')
xlim([1970 2100])
stamp('')
print('-dpng',[ppath 'Hist_SSP_ts_inputs_raw.png'])

%% diff from 1990 sep 
figure(5)
subplot(2,2,1)
plot(yH,GHist_Tp-GHistTp90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GSSP126_Tp-GHistTp90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_Tp-GHistTp90,'r','LineWidth',1); hold on;
title('Tp')

subplot(2,2,2)
plot(yH,GHist_Tb-GHistTb90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GSSP126_Tb-GHistTb90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_Tb-GHistTb90,'r','LineWidth',1); hold on;
title('Tb')

subplot(2,2,3)
plot(yH,GHist_mz-GHistMz90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GSSP126_mz-GHistMz90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_mz-GHistMz90,'r','LineWidth',1); hold on;
ylabel('Biomass (g m^-^2) difference from 1990-2000')
title('Mesoz')

subplot(2,2,4)
plot(yH,GHist_det-GHistDet90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GSSP126_det-GHistDet90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_det-GHistDet90,'r','LineWidth',1); hold on;
title('Det')
ylabel('Flux (g m^-^2 s^-^1) difference from 1990-2000')
stamp('gfdl ')
print('-dpng',[ppath 'gfdl_Hist_SSP_ts_inputs_diff_1990_2000.png'])


%IPSL
figure(6)
subplot(2,2,1)
plot(yH,IHist_Tp-IHistTp90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,ISSP126_Tp-IHistTp90,'b','LineWidth',1); hold on;
plot(yS,ISSP585_Tp-IHistTp90,'r','LineWidth',1); hold on;
title('Tp')

subplot(2,2,2)
plot(yH,IHist_Tb-IHistTb90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,ISSP126_Tb-IHistTb90,'b','LineWidth',1); hold on;
plot(yS,ISSP585_Tb-IHistTb90,'r','LineWidth',1); hold on;
title('Tb')

subplot(2,2,3)
plot(yH,IHist_mz-IHistMz90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,ISSP126_mz-IHistMz90,'b','LineWidth',1); hold on;
plot(yS,ISSP585_mz-IHistMz90,'r','LineWidth',1); hold on;
ylabel('Biomass (g m^-^2) difference from 1990-2000')
title('Mesoz')

subplot(2,2,4)
plot(yH,IHist_det-IHistDet90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,ISSP126_det-IHistDet90,'b','LineWidth',1); hold on;
plot(yS,ISSP585_det-IHistDet90,'r','LineWidth',1); hold on;
title('Det')
ylabel('Flux (g m^-^2 s^-^1) difference from 1990-2000')
stamp('ipsl ')
print('-dpng',[ppath 'ipsl_Hist_SSP_ts_inputs_diff_1990_2000.png'])

%% diff from 1990 together
figure(3)
subplot(2,2,1)
plot(yH,GHist_Tp-GHistTp90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GSSP126_Tp-GHistTp90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_Tp-GHistTp90,'r','LineWidth',1); hold on;
plot(yH,IHist_Tp-IHistTp90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,ISSP126_Tp-IHistTp90,'c','LineWidth',1); hold on;
plot(yS,ISSP585_Tp-IHistTp90,'m','LineWidth',1); hold on;
title('Tp')
legend('Ghist','G126','G585','Ihist','I126','I585')
legend('location','northwest')

subplot(2,2,2)
plot(yH,GHist_Tb-GHistTb90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GSSP126_Tb-GHistTb90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_Tb-GHistTb90,'r','LineWidth',1); hold on;
plot(yH,IHist_Tb-IHistTb90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,ISSP126_Tb-IHistTb90,'c','LineWidth',1); hold on;
plot(yS,ISSP585_Tb-IHistTb90,'m','LineWidth',1); hold on;
title('Tb')

subplot(2,2,3)
plot(yH,GHist_mz-GHistMz90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GSSP126_mz-GHistMz90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_mz-GHistMz90,'r','LineWidth',1); hold on;
plot(yH,IHist_mz-IHistMz90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,ISSP126_mz-IHistMz90,'c','LineWidth',1); hold on;
plot(yS,ISSP585_mz-IHistMz90,'m','LineWidth',1); hold on;
ylabel('Biomass (g m^)-^2) difference from 1990-2000')
title('Mesoz')

subplot(2,2,4)
plot(yH,GHist_det-GHistDet90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GSSP126_det-GHistDet90,'b','LineWidth',1); hold on;
plot(yS,GSSP585_det-GHistDet90,'r','LineWidth',1); hold on;
plot(yH,IHist_det-IHistDet90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,ISSP126_det-IHistDet90,'c','LineWidth',1); hold on;
plot(yS,ISSP585_det-IHistDet90,'m','LineWidth',1); hold on;
title('Det')
ylabel('Flux (g m^-^2 s^-^1) difference from 1990-2000')
stamp('')
print('-dpng',[ppath 'Hist_SSP_ts_inputs_diff_1990_2000.png'])

