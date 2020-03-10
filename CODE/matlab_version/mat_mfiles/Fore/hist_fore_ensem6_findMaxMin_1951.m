% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Find ensemble end members with max and min change relative to their 1951
% value

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
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat']);

%% Ensemble parameter sets
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);

lpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];

load([lpath 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext');

%% ts
%In original saved file
% y1 = 1860+(1/12):(1/12):2005;
% y2 = 2005+(1/12):(1/12):2100;
% y = [y1 y2];

% SOMETHING WRONG WITH MOVING MEAN OF #28 AROUND 1949
% FIX TS OF SMALL FISH AT T=1081
%tF[1081,28]=1.94e+35
%tP[1081,28]=1.94e+35
%tD[1081,28]=1.94e+35
%tA[1081,28]=5.83e+35
hTsF(28,1081) = (hTsF(28,1080) + hTsF(28,1082))/2;
hTsP(28,1081) = (hTsP(28,1080) + hTsP(28,1082))/2;
hTsD(28,1081) = (hTsD(28,1080) + hTsD(28,1082))/2;

HF = hTsF + hTmF;
HP = hTsP + hTmP + hTlP;
HD = hTsD + hTmD + hTlD;
HA = HF + HP + HD;

FF = fTsF + fTmF;
FP = fTsP + fTmP + fTlP;
FD = fTsD + fTmD + fTlD;
FA = FF + FP + FD;

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];
tB = [hTB fTB];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];

mtF = mean(tF);
mtP = mean(tP);
mtD = mean(tD);
mtA = mean(tA);
mtB = mean(tB);

%% Variability 
hyr = find(y>1950 & y<=2000);
fyr = find(y>2050 & y<=2100);

vstats(1,1) = mean(var(tA(:,hyr),0,2));
vstats(1,2) = mean(var(tA(:,fyr),0,2));
vstats(2,1) = mean(var(tF(:,hyr),0,2));
vstats(2,2) = mean(var(tF(:,fyr),0,2));
vstats(3,1) = mean(var(tP(:,hyr),0,2));
vstats(3,2) = mean(var(tP(:,fyr),0,2));
vstats(4,1) = mean(var(tD(:,hyr),0,2));
vstats(4,2) = mean(var(tD(:,fyr),0,2));
vstats(5,1) = mean(var(tB(:,hyr),0,2));
vstats(5,2) = mean(var(tB(:,fyr),0,2));

vstats(1,3) = var(mtA(:,hyr));
vstats(1,4) = var(mtA(:,fyr));
vstats(2,3) = var(mtF(:,hyr));
vstats(2,4) = var(mtF(:,fyr));
vstats(3,3) = var(mtP(:,hyr));
vstats(3,4) = var(mtP(:,fyr));
vstats(4,3) = var(mtD(:,hyr));
vstats(4,4) = var(mtD(:,fyr));
vstats(5,3) = var(mtB(:,hyr));
vstats(5,4) = var(mtB(:,fyr));

Stab = array2table(vstats,'VariableNames',{'HvarAll','FvarAll','HvarMean','FvarMean'},...
    'RowNames',{'All','F','P','D','B'});
writetable(Stab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_Var.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Stab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_Var.csv'],...
    'Delimiter',',','WriteRowNames',true)

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
pstats(9,1)  = find(diffB==max(diffB));
pstats(10,1) = find(diffB==min(diffB));

pstats(1,2:7) = red_params(pstats(1,1),:);
pstats(2,2:7) = red_params(pstats(2,1),:);
pstats(3,2:7) = red_params(pstats(3,1),:);
pstats(4,2:7) = red_params(pstats(4,1),:);
pstats(5,2:7) = red_params(pstats(5,1),:);
pstats(6,2:7) = red_params(pstats(6,1),:);
pstats(7,2:7) = red_params(pstats(7,1),:);
pstats(8,2:7) = red_params(pstats(8,1),:);
pstats(9,2:7)  = red_params(pstats(9,1),:);
pstats(10,2:7) = red_params(pstats(10,1),:);

Ttab = array2table(tstats,'VariableNames',{'max','mean','min'},...
    'RowNames',{'All','F','P','D','B'});

Ptab = array2table(pstats,'VariableNames',['sim',ptext],...
    'RowNames',{'maxAll','minAll','maxF','minF','maxP','minP','maxD','minD',...
    'maxB','minB'});

writetable(Ttab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffs.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ttab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffs.csv'],...
    'Delimiter',',','WriteRowNames',true)

writetable(Ptab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffSims.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ptab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MaxMinDiffSims.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% save
save([epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims.mat'],...
    'Stab','vstats','Ttab','tstats','Ptab','pstats');
save([lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims.mat'],...
    'Stab','vstats','Ttab','tstats','Ptab','pstats');

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
ylabel('Change in biomass (g m^-^2) relative to 1951')
title('All')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinAll.png'])

%% F
figure(2)
plot(y,dtF(pstats(3,1),:),'r','LineWidth',2); hold on;
plot(y,dmF,'k','LineWidth',2); hold on;
plot(y,dtF(pstats(4,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in biomass (g m^-^2) relative to 1951')
title('Forage fish')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinF.png'])

%% P
figure(3)
plot(y,dtP(pstats(5,1),:),'r','LineWidth',2); hold on;
plot(y,dmP,'k','LineWidth',2); hold on;
plot(y,dtP(pstats(6,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in biomass (g m^-^2) relative to 1951')
title('Large pelagic fish')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinP.png'])

%% D
figure(4)
plot(y,dtD(pstats(7,1),:),'r','LineWidth',2); hold on;
plot(y,dmD,'k','LineWidth',2); hold on;
plot(y,dtD(pstats(8,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in biomass (g m^-^2) relative to 1951')
title('Demersal fish')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinD.png'])

%% B
figure(10)
plot(y,dtB(pstats(9,1),:),'r','LineWidth',2); hold on;
plot(y,dmB,'k','LineWidth',2); hold on;
plot(y,dtB(pstats(10,1),:),'b','LineWidth',2); hold on;
xlim([1950 2100])
%ylim([0.425 0.725])
legend('max','mean','min')
legend('location','southwest')
ylabel('Change in biomass (g m^-^2) relative to 1951')
title('Benthos')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_MaxMinB.png'])

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
xspec = [39;33;35;36;26];
xtext = {'max \DeltaAll','max \DeltaF','max \DeltaP','max \DeltaD','max \DeltaB'};

figure(5)
subplot(3,2,5)
plot(y,dtA(xspec,:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('All')

subplot(3,2,1)
plot(y,dtF(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.3 0.3])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(3,2,2)
plot(y,dtP(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.7 0.1])
title('Large pelagics')

subplot(3,2,3)
plot(y,dtD(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.25 0.05])
title('Demersals')
legend(xtext)
legend('location','southwest')

subplot(3,2,4)
plot(y,dtB(xspec,:)); hold on;
xlim([1950 2100])
ylim([-0.2 0.2])
title('Benthos')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_max_ts.png'])


%% Min
nspec = [2;35;3;17;33];
ntext = {'min \DeltaAll','min \DeltaF','min \DeltaP','min \DeltaD','min \DeltaB'};

figure(6)
subplot(3,2,5)
plot(y,dtA(nspec,:)); hold on;
xlim([1950 2100])
ylim([-1 0.3])
title('All')

subplot(3,2,1)
plot(y,dtF(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.3 0.3])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(3,2,2)
plot(y,dtP(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.7 0.1])
legend(ntext)
legend('location','southwest')
title('Large pelagics')

subplot(3,2,3)
plot(y,dtD(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.25 0.05])
title('Demersals')

subplot(3,2,4)
plot(y,dtB(nspec,:)); hold on;
xlim([1950 2100])
ylim([-0.2 0.2])
title('Benthos')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_min_ts.png'])

%% all types together
figure(7)
subplot(3,2,5)
plot(y,dtA(xspec(1),:)); hold on;
plot(y,dtF(xspec(1),:)); hold on;
plot(y,dtP(xspec(1),:)); hold on;
plot(y,dtD(xspec(1),:)); hold on;
plot(y,dtB(xspec(1),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('max \DeltaAll')

subplot(3,2,1)
plot(y,dtA(xspec(2),:)); hold on;
plot(y,dtF(xspec(2),:)); hold on;
plot(y,dtP(xspec(2),:)); hold on;
plot(y,dtD(xspec(2),:)); hold on;
plot(y,dtB(xspec(2),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('max \DeltaF')

subplot(3,2,2)
plot(y,dtA(xspec(3),:)); hold on;
plot(y,dtF(xspec(3),:)); hold on;
plot(y,dtP(xspec(3),:)); hold on;
plot(y,dtD(xspec(3),:)); hold on;
plot(y,dtB(xspec(3),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('max \DeltaP')

subplot(3,2,3)
plot(y,dtA(xspec(4),:)); hold on;
plot(y,dtF(xspec(4),:)); hold on;
plot(y,dtP(xspec(4),:)); hold on;
plot(y,dtD(xspec(4),:)); hold on;
plot(y,dtB(xspec(4),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('max \DeltaD')
xlim([1950 2100])

subplot(3,2,4)
plot(y,dtA(xspec(5),:)); hold on;
plot(y,dtF(xspec(5),:)); hold on;
plot(y,dtP(xspec(5),:)); hold on;
plot(y,dtD(xspec(5),:)); hold on;
plot(y,dtB(xspec(5),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('max \DeltaB')
xlim([1950 2100])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_max_ts_v2.png'])

%%
figure(8)
subplot(3,2,5)
plot(y,dtA(nspec(1),:)); hold on;
plot(y,dtF(nspec(1),:)); hold on;
plot(y,dtP(nspec(1),:)); hold on;
plot(y,dtD(nspec(1),:)); hold on;
plot(y,dtB(nspec(1),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('min \DeltaAll')

subplot(3,2,1)
plot(y,dtA(nspec(2),:)); hold on;
plot(y,dtF(nspec(2),:)); hold on;
plot(y,dtP(nspec(2),:)); hold on;
plot(y,dtD(nspec(2),:)); hold on;
plot(y,dtB(nspec(2),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('min \DeltaF')

subplot(3,2,2)
plot(y,dtA(nspec(3),:)); hold on;
plot(y,dtF(nspec(3),:)); hold on;
plot(y,dtP(nspec(3),:)); hold on;
plot(y,dtD(nspec(3),:)); hold on;
plot(y,dtB(nspec(3),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('min \DeltaP')

subplot(3,2,3)
plot(y,dtA(nspec(4),:)); hold on;
plot(y,dtF(nspec(4),:)); hold on;
plot(y,dtP(nspec(4),:)); hold on;
plot(y,dtD(nspec(4),:)); hold on;
plot(y,dtB(nspec(4),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('min \DeltaD')
xlim([1950 2100])

subplot(3,2,4)
plot(y,dtA(nspec(5),:)); hold on;
plot(y,dtF(nspec(5),:)); hold on;
plot(y,dtP(nspec(5),:)); hold on;
plot(y,dtD(nspec(5),:)); hold on;
plot(y,dtB(nspec(5),:)); hold on;
xlim([1950 2100])
ylim([-1 0.4])
title('min \DeltaB')
xlim([1950 2100])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_min_ts_v2.png'])





