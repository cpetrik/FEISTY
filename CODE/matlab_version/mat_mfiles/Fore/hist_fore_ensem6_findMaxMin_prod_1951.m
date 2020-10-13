% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Find ensemble end members with max and min change in production relative 
% to their 1951 value

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Original parameters
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
load([fpath 'ESM2M_Hist_Fore/Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
    'HF_tamean','HP_tamean','HD_tamean','HB_tamean',...
    'FF_tamean','FP_tamean','FD_tamean','FB_tamean',...
    'HA_tamean','FA_tamean');

%% Ensemble parameter sets
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'haTsF','haTsP','haTsD','haTmF','haTmP','haTmD','haTB','haTlP','haTlD');
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'faTsF','faTsP','faTsD','faTmF','faTmP','faTmD','faTB','faTlP','faTlD');

%aT=tamean=nansum(prod.*area,1);

% lpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',...
%     'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext');

%% calculate: 
%NO benthic production = benthic biomass * detritus flux * benthic efficiency
%benthic production = detritus flux * benthic efficiency
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([bpath 'cobalt_det_ts.mat'],'det_mean_hist','det_mean_fore',...
    'mo_hist','mo_fore');
det_hist = det_mean_hist' * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore' * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

% HB = haTB .* repmat(det_hist,43,1) .* 0.075;
% FB = faTB .* repmat(det_fore,43,1) .* 0.075;
% HB_tam = HB_tamean .* det_hist .* 0.075;
% FB_tam = FB_tamean .* det_fore .* 0.075;

HB = repmat(det_hist,43,1) .* 0.075;
FB = repmat(det_fore,43,1) .* 0.075;
HB_tam = det_hist .* 0.075;
FB_tam = det_fore .* 0.075;

%% ts
%In original saved file
y1 = 1860+(1/12):(1/12):2005;
y2 = 2005+(1/12):(1/12):2100;
y = [y1 y2];

HF = haTsF + haTmF;
HP = haTsP + haTmP + haTlP;
HD = haTsD + haTmD + haTlD;
HA = HF + HP + HD;

FF = faTsF + faTmF;
FP = faTsP + faTmP + faTlP;
FD = faTsD + faTmD + faTlD;
FA = FF + FP + FD;

tForig = [HF_tamean FF_tamean];
tPorig = [HP_tamean FP_tamean];
tDorig = [HD_tamean FD_tamean];
tAorig = [HA_tamean FA_tamean];
tBorig = [HB_tam FB_tam];

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];
tB = [HB FB];
%tB = [haTB faTB];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];
tB = [tB; tBorig];

mtF = mean(tF);
mtP = mean(tP);
mtD = mean(tD);
mtA = mean(tA);
mtB = mean(tB);

%% Variability - this means nothing without removing the trend first
hyr = find(y>1950 & y<=2000);
fyr = find(y>2050 & y<=2100);

% vstats(1,1) = mean(var(tA(:,hyr),0,2));
% vstats(1,2) = mean(var(tA(:,fyr),0,2));
% vstats(2,1) = mean(var(tF(:,hyr),0,2));
% vstats(2,2) = mean(var(tF(:,fyr),0,2));
% vstats(3,1) = mean(var(tP(:,hyr),0,2));
% vstats(3,2) = mean(var(tP(:,fyr),0,2));
% vstats(4,1) = mean(var(tD(:,hyr),0,2));
% vstats(4,2) = mean(var(tD(:,fyr),0,2));
% vstats(5,1) = mean(var(tB(:,hyr),0,2));
% vstats(5,2) = mean(var(tB(:,fyr),0,2));
% 
% vstats(1,3) = var(mtA(:,hyr));
% vstats(1,4) = var(mtA(:,fyr));
% vstats(2,3) = var(mtF(:,hyr));
% vstats(2,4) = var(mtF(:,fyr));
% vstats(3,3) = var(mtP(:,hyr));
% vstats(3,4) = var(mtP(:,fyr));
% vstats(4,3) = var(mtD(:,hyr));
% vstats(4,4) = var(mtD(:,fyr));
% vstats(5,3) = var(mtB(:,hyr));
% vstats(5,4) = var(mtB(:,fyr));
% 
% Stab = array2table(vstats,'VariableNames',{'HvarAll','FvarAll','HvarMean','FvarMean'},...
%     'RowNames',{'All','F','P','D','B'});
% writetable(Stab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_Var_prod.csv'],...
%     'Delimiter',',','WriteRowNames',true)
% writetable(Stab,[dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_Var_prod.csv'],...
%     'Delimiter',',','WriteRowNames',true)

%% Difference
diffA = tA(:,hyr(1)) - tA(:,end);
diffF = tF(:,hyr(1)) - tF(:,end);
diffP = tP(:,hyr(1)) - tP(:,end);
diffD = tD(:,hyr(1)) - tD(:,end);
diffB = tB(:,hyr(1)) - tB(:,end);

tstats(1,1) = max(diffA);
tstats(1,2) = mean(diffA);
tstats(1,3) = min(diffA);
tstats(2,1) = max(diffF);
tstats(2,2) = mean(diffF);
tstats(2,3) = min(diffF);
tstats(3,1) = max(diffP);
tstats(3,2) = mean(diffP);
tstats(3,3) = min(diffP);
tstats(4,1) = max(diffD);
tstats(4,2) = mean(diffD);
tstats(4,3) = min(diffD);
tstats(5,1) = max(diffB);
tstats(5,2) = mean(diffB);
tstats(5,3) = min(diffB);

pstats(1,1) = find(diffA==max(diffA));
pstats(2,1) = find(diffA==min(diffA));
pstats(3,1) = find(diffF==max(diffF));
pstats(4,1) = find(diffF==min(diffF));
pstats(5,1) = find(diffP==max(diffP));
pstats(6,1) = find(diffP==min(diffP));
pstats(7,1) = find(diffD==max(diffD));
pstats(8,1) = find(diffD==min(diffD));
% pstats(9,1)  = find(diffB==max(diffB));
% pstats(10,1) = find(diffB==min(diffB));

pstats(1,2:7) = red_params(pstats(1,1),:);
pstats(2,2:7) = red_params(pstats(2,1),:);
pstats(3,2:7) = red_params(pstats(3,1),:);
pstats(4,2:7) = red_params(pstats(4,1),:);
pstats(5,2:7) = red_params(pstats(5,1),:);
pstats(6,2:7) = red_params(pstats(6,1),:);
pstats(7,2:7) = red_params(pstats(7,1),:);
pstats(8,2:7) = red_params(pstats(8,1),:);
% pstats(9,2:7)  = red_params(pstats(9,1),:);
% pstats(10,2:7) = red_params(pstats(10,1),:);

Ttab = array2table(tstats,'VariableNames',{'max','mean','min'},...
    'RowNames',{'All','F','P','D','B'});

Ptab = array2table(pstats,'VariableNames',['sim',ptext],...
    'RowNames',{'maxAll','minAll','maxF','minF','maxP','minP','maxD','minD'});%,'maxB','minB'});

writetable(Ttab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffs_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ttab,[dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffs_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)

writetable(Ptab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffSims_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ptab,[dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffSims_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% save
save([epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],...
    'Stab','vstats','Ttab','tstats','Ptab','pstats');
save([dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],...
    'Stab','vstats','Ttab','tstats','Ptab','pstats');

%% Percent Difference
% percent difference from 1951?
pdF = (tF(:,hyr(1)) - tF(:,end)) ./ tF(:,hyr(1)); 
pdP = (tP(:,hyr(1)) - tP(:,end)) ./ tP(:,hyr(1));
pdD = (tD(:,hyr(1)) - tD(:,end)) ./ tD(:,hyr(1)); 
pdA = (tA(:,hyr(1)) - tA(:,end)) ./ tA(:,hyr(1));
pdB = (tB(:,hyr(1)) - tB(:,end)) ./ tB(:,hyr(1));

mmstat(1,1) = max(pdA);
mmstat(1,2) = mean(pdA);
mmstat(1,3) = min(pdA);
mmstat(2,1) = max(pdF);
mmstat(2,2) = mean(pdF);
mmstat(2,3) = min(pdF);
mmstat(3,1) = max(pdP);
mmstat(3,2) = mean(pdP);
mmstat(3,3) = min(pdP);
mmstat(4,1) = max(pdD);
mmstat(4,2) = mean(pdD);
mmstat(4,3) = min(pdD);
mmstat(5,1) = max(pdB);
mmstat(5,2) = mean(pdB);
mmstat(5,3) = min(pdB);

pdstat(1,1) = find(pdA==max(pdA));
pdstat(2,1) = find(pdA==min(pdA));
pdstat(3,1) = find(pdF==max(pdF));
pdstat(4,1) = find(pdF==min(pdF));
pdstat(5,1) = find(pdP==max(pdP));
pdstat(6,1) = find(pdP==min(pdP));
pdstat(7,1) = find(pdD==max(pdD));
pdstat(8,1) = find(pdD==min(pdD));
% pdstat(9,1)  = find(pdB==max(pdB));
% pdstat(10,1) = find(pdB==min(pdB));

pdstat(1,2:7) = red_params(pdstat(1,1),:);
pdstat(2,2:7) = red_params(pdstat(2,1),:);
pdstat(3,2:7) = red_params(pdstat(3,1),:);
pdstat(4,2:7) = red_params(pdstat(4,1),:);
pdstat(5,2:7) = red_params(pdstat(5,1),:);
pdstat(6,2:7) = red_params(pdstat(6,1),:);
pdstat(7,2:7) = red_params(pdstat(7,1),:);
pdstat(8,2:7) = red_params(pdstat(8,1),:);
% pdstat(9,2:7)  = red_params(pdstat(9,1),:);
% pdstat(10,2:7) = red_params(pdstat(10,1),:);

Mtab = array2table(mmstat,'VariableNames',{'max','mean','min'},...
    'RowNames',{'All','F','P','D','B'});

Dtab = array2table(pdstat,'VariableNames',['sim',ptext],...
    'RowNames',{'maxAll','minAll','maxF','minF','maxP','minP','maxD','minD',...
    });%'maxB','minB'});

writetable(Mtab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinpPdiffs_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Mtab,[dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinPdiffs_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)

writetable(Dtab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinPdiffSims_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Dtab,[dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinPdiffSims_prod.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% save
save([epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],...
    'Mtab','mmstat','Dtab','pdstat','-append');
save([dpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],...
    'Mtab','mmstat','Dtab','pdstat','-append');

%% Line color order
cm21=[1 0 0;...    %r
    0 0 0.75;...  %b
    0 0.7 0;...   %g
    0 0 0;...     %black
    0.5 0 0;...   %maroon
    1 0 1;...     %m
    0.5 0 1;...   %purple
    0 0.5 0.75;...   %med blue
    0/255 206/255 209/255;... %turq
    0 1 1;...     %c
    0 0.5 0;...   %dk green
    0.5 0.5 0;... %tan/army
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255;... %peach
    1 0.5 0];  %orange
    
set(groot,'defaultAxesColorOrder',cm21);

%% moving means
mmtF = movmean(tF,12,2); %[1081,28]=1.94e+35
mmtP = movmean(tP,12,2); %[1081,28]=1.94e+35
mmtD = movmean(tD,12,2); %[1081,28]=1.94e+35
mmtA = movmean(tA,12,2); %[1081,28]=5.83e+35
mmtB = movmean(tB,12,2);

mmmF = movmean(mtF,12);
mmmP = movmean(mtP,12);
mmmD = movmean(mtD,12);
mmmA = movmean(mtA,12);
mmmB = movmean(mtB,12);

dtF = mmtF - mmtF(:,hyr(1));
dtP = mmtP - mmtP(:,hyr(1));
dtD = mmtD - mmtD(:,hyr(1));
dtA = mmtA - mmtA(:,hyr(1));
dtB = mmtB - mmtB(:,hyr(1));

dmF = mmmF - mmmF(:,hyr(1));
dmP = mmmP - mmmP(:,hyr(1));
dmD = mmmD - mmmD(:,hyr(1));
dmA = mmmA - mmmA(:,hyr(1));
dmB = mmmB - mmmB(:,hyr(1));

%% Individual Plots 
% All
figure(1)
%subplot(3,3,1)
plot(y,dtA(pstats(1,1),:),'r','LineWidth',2); hold on;
plot(y,dmA,'k','LineWidth',2); hold on;
plot(y,dtA(pstats(2,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in production (g d^-^1) relative to 1951')
title('All')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinAll_prod.png'])

%% F
figure(2)
plot(y,dtF(pstats(3,1),:),'r','LineWidth',2); hold on;
plot(y,dmF,'k','LineWidth',2); hold on;
plot(y,dtF(pstats(4,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in production (g d^-^1) relative to 1951')
title('Forage fish')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinF_prod.png'])

%% P
figure(3)
plot(y,dtP(pstats(5,1),:),'r','LineWidth',2); hold on;
plot(y,dmP,'k','LineWidth',2); hold on;
plot(y,dtP(pstats(6,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in production (g d^-^1) relative to 1951')
title('Large pelagic fish')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinP_prod.png'])

%% D
figure(4)
plot(y,dtD(pstats(7,1),:),'r','LineWidth',2); hold on;
plot(y,dmD,'k','LineWidth',2); hold on;
plot(y,dtD(pstats(8,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in production (g d^-^1) relative to 1951')
title('Demersal fish')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinD_prod.png'])

%% B
% figure(10)
% plot(y,dtB(pstats(9,1),:),'r','LineWidth',2); hold on;
% plot(y,dmB,'k','LineWidth',2); hold on;
% plot(y,dtB(pstats(10,1),:),'b','LineWidth',2); hold on;
% xlim([1950 2100])
% %ylim([0.425 0.725])
% legend('max','mean','min')
% legend('location','southwest')
% ylabel('Change in production (g d^-^1) relative to 1951')
% title('Benthos')
% print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinB_prod.png'])

%%
% Line color order
cm5=[0.5 0 1;...   %purple
    1 0 0;...     %r
    0 0 0.75;...  %b
    0 0.5 0;... %dk green
    188/255 143/255 143/255; ... %rosy brown
    0 0 0]; 

set(groot,'defaultAxesColorOrder',cm5);

%% Max
xspec = pstats(1:2:end,1);
xtext = {'max \DeltaAll','max \DeltaF','max \DeltaP','max \DeltaD','max \DeltaB'};

figure(5)
subplot(3,2,4)
plot(y,dtA(xspec,:)); hold on;
xlim([1950 2100])
ylim([-7e-3 3e-3])
title('All')

subplot(3,2,1)
plot(y,dtF(xspec,:)); hold on;
xlim([1950 2100])
ylim([-3e-3 2e-3])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(3,2,2)
plot(y,dtP(xspec,:)); hold on;
xlim([1950 2100])
ylim([-4e-3 1e-3])
title('Large pelagics')

subplot(3,2,3)
plot(y,dtD(xspec,:)); hold on;
xlim([1950 2100])
ylim([-5.5e-4 1e-4])
title('Demersals')
legend(xtext)
legend('location','southwest')

subplot(3,2,5)
plot(y,dtB(xspec,:)); hold on;
xlim([1950 2100])
ylim([-6e-4 1e-4])
title('Benthos')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_max_ts_prod.png'])


%% Min
nspec = pstats(2:2:end,1);
ntext = {'min \DeltaAll','min \DeltaF','min \DeltaP','min \DeltaD','min \DeltaB'};

figure(6)
subplot(3,2,4)
plot(y,dtA(nspec,:)); hold on;
xlim([1950 2100])
ylim([-6e-3 2e-3])
title('All')

subplot(3,2,1)
plot(y,dtF(nspec,:)); hold on;
xlim([1950 2100])
ylim([-3e-3 2e-3])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(3,2,2)
plot(y,dtP(nspec,:)); hold on;
xlim([1950 2100])
ylim([-1.5e-3 0.5e-3])
title('Large pelagics')

subplot(3,2,3)
plot(y,dtD(nspec,:)); hold on;
xlim([1950 2100])
ylim([-5e-4 1e-4])
legend(ntext)
legend('location','southwest')
title('Demersals')

subplot(3,2,5)
plot(y,dtB(nspec,:)); hold on;
xlim([1950 2100])
ylim([-5e-4 2e-4])
title('Benthos')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_min_ts_prod.png'])

%% ------------------- P Diffs ------------------------------------------
% moving means
% mmtF = movmean(tF,12,2); 
% mmtP = movmean(tP,12,2); 
% mmtD = movmean(tD,12,2); 
% mmtA = movmean(tA,12,2); 
% mmtB = movmean(tB,12,2);

%dtF = mmtF - mmtF(:,hyr(1));
%percent difference from 1951
mpdF = (-mmtF(:,hyr(1)) + mmtF) ./ mmtF(:,hyr(1)); 
mpdP = (-mmtP(:,hyr(1)) + mmtP) ./ mmtP(:,hyr(1));
mpdD = (-mmtD(:,hyr(1)) + mmtD) ./ mmtD(:,hyr(1)); 
mpdA = (-mmtA(:,hyr(1)) + mmtA) ./ mmtA(:,hyr(1));
mpdB = (-mmtB(:,hyr(1)) + mmtB) ./ mmtB(:,hyr(1));


%% Max
xspec = pdstat(1:2:end,1);
xtext = {'max %\DeltaAll','max %\DeltaF','max %\DeltaP','max %\DeltaD',...
    'max %\DeltaB'};

figure(11)
subplot(3,2,4)
plot(y,mpdA(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.2 0.1])
title('All')

subplot(3,2,1)
plot(y,mpdF(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.2 0.1])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(3,2,2)
plot(y,mpdP(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.4 0.1])
title('Large pelagics')

subplot(3,2,3)
plot(y,mpdD(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.25 0.05])
title('Demersals')
% legend(xtext)
% legend('location','southwest')

subplot(3,2,5)
plot(y,mpdB(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.4 0.2])
title('Benthos')

subplot(3,2,6)
plot(y,mpdD(xspec,:)); hold on;
xlim([1950 2100])
ylim([1 2])
legend(xtext)
legend('location','southwest')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_max_ts_prod_pdiff.png'])


%% Min
nspec = pdstat(2:2:end,1);
ntext = {'min %\DeltaAll','min %\DeltaF','min %\DeltaP','min %\DeltaD',...
    'min %\DeltaB'};

figure(12)
subplot(3,2,4)
plot(y,mpdA(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.2 0.1])
title('All')

subplot(3,2,1)
plot(y,mpdF(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.15 0.1])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(3,2,2)
plot(y,mpdP(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.35 0.1])
title('Large pelagics')

subplot(3,2,3)
plot(y,mpdD(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.25 0.05])
% legend(ntext)
% legend('location','southwest')
title('Demersals')

subplot(3,2,5)
plot(y,mpdB(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.4 0.2])
title('Benthos')

subplot(3,2,6)
plot(y,mpdD(xspec,:)); hold on;
xlim([1950 2100])
ylim([1 2])
legend(ntext)
legend('location','southwest')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_min_ts_prod_pdiff.png'])



