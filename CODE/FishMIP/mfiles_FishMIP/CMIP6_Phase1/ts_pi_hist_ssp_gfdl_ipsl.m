% Visualize time series output of FEISTY forced by CMIP6

clear all
close all

%% Fish data
cfile = 'Dc_Lam579_enc70-b200_m440-b175-k086_c20-b250_D080_A050_nmort1_BE10_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%gfdl
gpath=['/Volumes/FEISTY/NC/FishMIP/GFDL_CMIP6/' cfile '/'];
% load([gpath 'Means_PreIndust_' cfile '.mat'],...
%     'GPreFts','GPrePts','GPreDts','GPreBts');
load([gpath 'Means_Hist_' cfile '.mat'],...
    'GHistFts','GHistPts','GHistDts','GHistBts');
load([gpath 'Means_SSP126_2090-2100_' cfile '.mat'],...
    'GS126Fts','GS126Pts','GS126Dts','GS126Bts');
load([gpath 'Means_SSP585_2090-2100_' cfile '.mat'],...
    'GS585Fts','GS585Pts','GS585Dts','GS585Bts');

%ipsl
ipath=['/Volumes/FEISTY/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
% load([ipath 'Means_PreIndust_' cfile '.mat'],...
%     'IPreFts','IPrePts','IPreDts','IPreBts');
load([ipath 'Means_Hist_' cfile '.mat'],...
    'IHistFts','IHistPts','IHistDts','IHistBts');
load([ipath 'Means_SSP126_2090-2100_' cfile '.mat'],...
    'IS126Fts','IS126Pts','IS126Dts','IS126Bts');
load([ipath 'Means_SSP585_2090-2100_' cfile '.mat'],...
    'IS585Fts','IS585Pts','IS585Dts','IS585Bts');

%all fish
GPreAts  = GPreFts +GPrePts +GPreDts;
GHistAts = GHistFts+GHistPts+GHistDts;
GS126Ats = GS126Fts+GS126Pts+GS126Dts;
GS585Ats = GS585Fts+GS585Pts+GS585Dts;
IPreAts  = IPreFts +IPrePts +IPreDts;
IHistAts = IHistFts+IHistPts+IHistDts;
IS126Ats = IS126Fts+IS126Pts+IS126Dts;
IS585Ats = IS585Fts+IS585Pts+IS585Dts;

%% Combine Hist & SSP ts
% tF = [HF FF];
% tP = [HP FP];
% tD = [HD FD];
% tA = [HA FA];
% 
% 
% %% moving means
% mmtF = movmean(tF,12,2); 
% mmtP = movmean(tP,12,2); 
% mmtD = movmean(tD,12,2); 
% mmtA = movmean(tA,12,2); 

%% time
yP = 1950+(1/12):(1/12):2100;
yH = 1949+(1/12):(1/12):2014;
yS = 2014+(1/12):(1/12):2100;

%% pure biomass
figure(1)
subplot(2,2,1)
plot(yH,GHistFts,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Fts,'b','LineWidth',2); hold on;
plot(yS,GS585Fts,'r','LineWidth',2); hold on;
plot(yH,IHistFts,'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Fts,'--b','LineWidth',2); hold on;
plot(yS,IS585Fts,'--r','LineWidth',2); hold on;
ylabel('Biomass (g m^-^2)')
title('F')

subplot(2,2,2)
plot(yH,GHistPts,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Pts,'b','LineWidth',2); hold on;
plot(yS,GS585Pts,'r','LineWidth',2); hold on;
plot(yH,IHistPts,'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Pts,'--b','LineWidth',2); hold on;
plot(yS,IS585Pts,'--r','LineWidth',2); hold on;
ylabel('Biomass (g m^-^2)')
title('P')

subplot(2,2,3)
plot(yH,GHistDts,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Dts,'b','LineWidth',2); hold on;
plot(yS,GS585Dts,'r','LineWidth',2); hold on;
plot(yH,IHistDts,'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Dts,'--b','LineWidth',2); hold on;
plot(yS,IS585Dts,'--r','LineWidth',2); hold on;
ylabel('Biomass (g m^-^2)')
title('D')

subplot(2,2,4)
plot(yH,GHistAts,'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Ats,'b','LineWidth',2); hold on;
plot(yS,GS585Ats,'r','LineWidth',2); hold on;
plot(yH,IHistAts,'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Ats,'--b','LineWidth',2); hold on;
plot(yS,IS585Ats,'--r','LineWidth',2); hold on;
ylabel('Biomass (g m^-^2)')
title('All fish')
stamp('')
print('-dpng',[ppath 'Pre_Hist_SSP_ts_all_types_biom.png'])

%% diff from 1950
figure(2)
subplot(2,2,1)
plot(yH,GHistFts-GHistFts(1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Fts-GHistFts(1),'b','LineWidth',2); hold on;
plot(yS,GS585Fts-GHistFts(1),'r','LineWidth',2); hold on;
plot(yH,IHistFts-IHistFts(1),'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Fts-IHistFts(1),'--b','LineWidth',2); hold on;
plot(yS,IS585Fts-IHistFts(1),'--r','LineWidth',2); hold on;
title('F')

subplot(2,2,2)
plot(yH,GHistPts-GHistPts(1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Pts-GHistPts(1),'b','LineWidth',2); hold on;
plot(yS,GS585Pts-GHistPts(1),'r','LineWidth',2); hold on;
plot(yH,IHistPts-IHistPts(1),'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Pts-IHistPts(1),'--b','LineWidth',2); hold on;
plot(yS,IS585Pts-IHistPts(1),'--r','LineWidth',2); hold on;
title('P')

subplot(2,2,3)
plot(yH,GHistDts-GHistDts(1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Dts-GHistDts(1),'b','LineWidth',2); hold on;
plot(yS,GS585Dts-GHistDts(1),'r','LineWidth',2); hold on;
plot(yH,IHistDts-IHistDts(1),'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Dts-IHistDts(1),'--b','LineWidth',2); hold on;
plot(yS,IS585Dts-IHistDts(1),'--r','LineWidth',2); hold on;
ylabel('Biomass (g m^-^2) difference from 1950')
title('D')

subplot(2,2,4)
plot(yH,GHistAts-GHistAts(1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,GS126Ats-GHistAts(1),'b','LineWidth',2); hold on;
plot(yS,GS585Ats-GHistAts(1),'r','LineWidth',2); hold on;
plot(yH,IHistAts-IHistAts(1),'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(yS,IS126Ats-IHistAts(1),'--b','LineWidth',2); hold on;
plot(yS,IS585Ats-IHistAts(1),'--r','LineWidth',2); hold on;
title('All fish')
stamp('')
print('-dpng',[ppath 'Pre_Hist_SSP_ts_all_types_biom_diff.png'])

