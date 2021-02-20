% Visualize time series output of FEISTY forced by CMIP6
% Fig 3B in ms
% tcb relative to 1990-2000

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%gfdl
gpath=['/Volumes/MIP/NC/FishMIP/GFDL_CMIP6/' cfile '/'];
load([gpath 'Means_PreIndust_empHP_' cfile '.mat'],...
    'GPreFts','GPrePts','GPreDts','GPreBts');
load([gpath 'Means_Hist_2000-2010_' cfile '.mat'],...
    'GHistFts','GHistPts','GHistDts','GHistBts');
load([gpath 'Means_SSP126_empHP_2090-2100_' cfile '.mat'],...
    'GS126Fts','GS126Pts','GS126Dts','GS126Bts');
load([gpath 'Means_SSP585_empHP_2090-2100_' cfile '.mat'],...
    'GS585Fts','GS585Pts','GS585Dts','GS585Bts');

%ipsl
ipath=['/Volumes/MIP/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
load([ipath 'Means_PreIndust_empHP_' cfile '.mat'],...
    'IPreFts','IPrePts','IPreDts','IPreBts');
load([ipath 'Means_Hist_empHP_2000-2010_' cfile '.mat'],...
    'IHistFts','IHistPts','IHistDts','IHistBts');
load([ipath 'Means_SSP126_empHP_2090-2100_' cfile '.mat'],...
    'IS126Fts','IS126Pts','IS126Dts','IS126Bts');
load([ipath 'Means_SSP585_empHP_2090-2100_' cfile '.mat'],...
    'IS585Fts','IS585Pts','IS585Dts','IS585Bts');

%all consumers
GPreAts  = GPreFts +GPrePts +GPreDts +GPreBts;
GHistAts = GHistFts+GHistPts+GHistDts+GHistBts;
GS126Ats = GS126Fts+GS126Pts+GS126Dts+GS126Bts;
GS585Ats = GS585Fts+GS585Pts+GS585Dts+GS585Bts;
IPreAts  = IPreFts +IPrePts +IPreDts +IPreBts;
IHistAts = IHistFts+IHistPts+IHistDts+IHistBts;
IS126Ats = IS126Fts+IS126Pts+IS126Dts+IS126Bts;
IS585Ats = IS585Fts+IS585Pts+IS585Dts+IS585Bts;

%% time
yP = 1950+(1/12):(1/12):2101;
yH = 1950+(1/12):(1/12):2015;
yS = 2015+(1/12):(1/12):2101;

tid = find(yH>1990 & yH<=2000);

GPreFm90 =  mean(GPreFts(tid));
IPreFm90 =  mean(IPreFts(tid));
GHistFm90 = mean(GHistFts(tid));
IHistFm90 = mean(IHistFts(tid));

GPrePm90 =  mean(GPrePts(tid));
IPrePm90 =  mean(IPrePts(tid));
GHistPm90 = mean(GHistPts(tid));
IHistPm90 = mean(IHistPts(tid));

GPreDm90 =  mean(GPreDts(tid));
IPreDm90 =  mean(IPreDts(tid));
GHistDm90 = mean(GHistDts(tid));
IHistDm90 = mean(IHistDts(tid));

GPreAm90 =  mean(GPreAts(tid));
IPreAm90 =  mean(IPreAts(tid));
GHistAm90 = mean(GHistAts(tid));
IHistAm90 = mean(IHistAts(tid));

%% percent diff from 1990-2000
figure(1)
subplot(2,2,1)
plot(yP,GPreFts ./GPreFm90,'k','LineWidth',1); hold on;
plot(yH,GHistFts./GHistFm90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GS126Fts./GHistFm90,'b','LineWidth',1); hold on;
plot(yS,GS585Fts./GHistFm90,'r','LineWidth',1); hold on;
plot(yP,IPreFts ./IPreFm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yH,IHistFts./IHistFm90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,IS126Fts./IHistFm90,'c','LineWidth',1); hold on;
plot(yS,IS585Fts./IHistFm90,'m','LineWidth',1); hold on;
title('F')
xlim([1970 2100])

subplot(2,2,2)
plot(yP,GPrePts ./GPrePm90,'k','LineWidth',1); hold on;
plot(yH,GHistPts./GHistPm90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GS126Pts./GHistPm90,'b','LineWidth',1); hold on;
plot(yS,GS585Pts./GHistPm90,'r','LineWidth',1); hold on;
plot(yP,IPrePts ./IPrePm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yH,IHistPts./IHistPm90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,IS126Pts./IHistPm90,'c','LineWidth',1); hold on;
plot(yS,IS585Pts./IHistPm90,'m','LineWidth',1); hold on;
title('P')
xlim([1970 2100])

subplot(2,2,3)
plot(yP,GPreDts ./GPreDm90,'k','LineWidth',1); hold on;
plot(yH,GHistDts./GHistDm90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GS126Dts./GHistDm90,'b','LineWidth',1); hold on;
plot(yS,GS585Dts./GHistDm90,'r','LineWidth',1); hold on;
plot(yP,IPreDts ./IPreDm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yH,IHistDts./IHistDm90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,IS126Dts./IHistDm90,'c','LineWidth',1); hold on;
plot(yS,IS585Dts./IHistDm90,'m','LineWidth',1); hold on;
ylabel('Biomass (g m^-^2) difference from 1990-2000')
title('D')
xlim([1970 2100])

subplot(2,2,4)
plot(yP,GPreAts ./GPreAm90,'k','LineWidth',1); hold on;
plot(yH,GHistAts./GHistAm90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,GS126Ats./GHistAm90,'b','LineWidth',1); hold on;
plot(yS,GS585Ats./GHistAm90,'r','LineWidth',1); hold on;
plot(yP,IPreAts ./IPreAm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yH,IHistAts./IHistAm90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,IS126Ats./IHistAm90,'c','LineWidth',1); hold on;
plot(yS,IS585Ats./IHistAm90,'m','LineWidth',1); hold on;
title('All consumers')
xlim([1970 2100])
stamp('')
print('-dpng',[ppath 'Pre_Hist_SSP_ts_empHP_all_types_biom_pdiff_1990_2000.png'])


%%
figure(2)
% plot(yP,GPreAts ./GPreAm90,'k','LineWidth',1); hold on;
% plot(yH,GHistAts./GHistAm90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
% plot(yS,GS126Ats./GHistAm90,'b','LineWidth',1); hold on;
% plot(yS,GS585Ats./GHistAm90,'r','LineWidth',1); hold on;
% plot(yP,IPreAts ./IPreAm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
% plot(yH,IHistAts./IHistAm90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
% plot(yS,IS126Ats./IHistAm90,'c','LineWidth',1); hold on;
% plot(yS,IS585Ats./IHistAm90,'m','LineWidth',1); hold on;
% title('All consumers')
% ylabel('Biomass (g m^-^2) difference from 1990-2000')
% xlim([1970 2100])
% ylim([0.6 1.2])
% stamp('')
% print('-dpng',[ppath 'tcb_Pre_Hist_SSP_ts_empHP_biom_rdiff_1990_2000.png'])

plot(yP,(GPreAts-GPreAm90) ./GPreAm90,'k','LineWidth',1); hold on;
plot(yH,(GHistAts-GHistAm90)./GHistAm90,'color',[0.25 0.25 0.25],'LineWidth',1); hold on;
plot(yS,(GS126Ats-GHistAm90)./GHistAm90,'b','LineWidth',1); hold on;
plot(yS,(GS585Ats-GHistAm90)./GHistAm90,'r','LineWidth',1); hold on;
plot(yP,(IPreAts-IPreAm90) ./IPreAm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yH,(IHistAts-IHistAm90)./IHistAm90,'color',[0.75 0.75 0.75],'LineWidth',1); hold on;
plot(yS,(IS126Ats-IHistAm90)./IHistAm90,'c','LineWidth',1); hold on;
plot(yS,(IS585Ats-IHistAm90)./IHistAm90,'m','LineWidth',1); hold on;
title('All consumers')
ylabel('Biomass (g m^-^2) difference from 1990-2000')
xlim([1970 2100])
ylim([-0.40 0.20])
set(gca,'YTick',[-0.4:0.05:0.2])
stamp('')
print('-dpng',[ppath 'tcb_Pre_Hist_SSP_ts_empHP_biom_pdiff_1990_2000.png'])

%% diff from 1990
figure(5)
subplot(2,2,1)
plot(yP,GPreFts -GPreFm90,'k','LineWidth',1); hold on;
plot(yH,GHistFts-GHistFm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GS126Fts-GHistFm90,'b','LineWidth',1); hold on;
plot(yS,GS585Fts-GHistFm90,'r','LineWidth',1); hold on;
title('F')

subplot(2,2,2)
plot(yP,GPrePts -GPrePm90,'k','LineWidth',1); hold on;
plot(yH,GHistPts-GHistPm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GS126Pts-GHistPm90,'b','LineWidth',1); hold on;
plot(yS,GS585Pts-GHistPm90,'r','LineWidth',1); hold on;
title('P')

subplot(2,2,3)
plot(yP,GPreDts -GPreDm90,'k','LineWidth',1); hold on;
plot(yH,GHistDts-GHistDm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GS126Dts-GHistDm90,'b','LineWidth',1); hold on;
plot(yS,GS585Dts-GHistDm90,'r','LineWidth',1); hold on;
ylabel('Biomass (g m^-^2) difference from 1990-2000')
title('D')

subplot(2,2,4)
plot(yP,GPreAts -GPreAm90,'k','LineWidth',1); hold on;
plot(yH,GHistAts-GHistAm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,GS126Ats-GHistAm90,'b','LineWidth',1); hold on;
plot(yS,GS585Ats-GHistAm90,'r','LineWidth',1); hold on;
title('All consumers')
stamp('gfdl ')
print('-dpng',[ppath 'gfdl_Pre_Hist_SSP_ts_empHP_all_types_biom_diff_1990_2000.png'])


%IPSL
figure(6)
subplot(2,2,1)
plot(yP,IPreFts -IPreFm90,'k','LineWidth',1); hold on;
plot(yH,IHistFts-IHistFm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,IS126Fts-IHistFm90,'b','LineWidth',1); hold on;
plot(yS,IS585Fts-IHistFm90,'r','LineWidth',1); hold on;
title('F')

subplot(2,2,2)
plot(yP,IPrePts -IPrePm90,'k','LineWidth',1); hold on;
plot(yH,IHistPts-IHistPm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,IS126Pts-IHistPm90,'b','LineWidth',1); hold on;
plot(yS,IS585Pts-IHistPm90,'r','LineWidth',1); hold on;
title('P')

subplot(2,2,3)
plot(yP,IPreDts -IPreDm90,'k','LineWidth',1); hold on;
plot(yH,IHistDts-IHistDm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,IS126Dts-IHistDm90,'b','LineWidth',1); hold on;
plot(yS,IS585Dts-IHistDm90,'r','LineWidth',1); hold on;
ylabel('Biomass (g m^-^2) difference from 1990-2000')
title('D')

subplot(2,2,4)
plot(yP,IPreAts -IPreAm90,'k','LineWidth',1); hold on;
plot(yH,IHistAts-IHistAm90,'color',[0.5 0.5 0.5],'LineWidth',1); hold on;
plot(yS,IS126Ats-IHistAm90,'b','LineWidth',1); hold on;
plot(yS,IS585Ats-IHistAm90,'r','LineWidth',1); hold on;
title('All consumers')
stamp('ipsl ')
print('-dpng',[ppath 'ipsl_Pre_Hist_SSP_ts_empHP_all_types_biom_diff_1990_2000.png'])

