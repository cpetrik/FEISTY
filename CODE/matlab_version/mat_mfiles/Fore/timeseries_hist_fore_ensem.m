% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' ...
    'param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/'];

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
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem5_mid5_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem5_mid5_bestAIC_multFup_multPneg.mat']);

%% ts
%In original saved file
% y1 = 1860+(1/12):(1/12):2005;
% y2 = 2005+(1/12):(1/12):2100;
% y = [y1 y2];

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

%% Calc variability at 1900, 2000, and 2100

%Table with most & least fish
tabML(1,1) = find(FA(:,1140)==max(FA(:,1140)));
tabML(1,2) = find(FA(:,1140)==min(FA(:,1140)));
tabML(2,1) = find(FF(:,1140)==max(FF(:,1140)));
tabML(2,2) = find(FF(:,1140)==min(FF(:,1140)));
tabML(3,1) = find(FP(:,1140)==max(FP(:,1140)));
tabML(3,2) = find(FP(:,1140)==min(FP(:,1140)));
tabML(4,1) = find(FD(:,1140)==max(FD(:,1140)));
tabML(4,2) = find(FD(:,1140)==min(FD(:,1140)));

Ftab = array2table(tabML,'VariableNames',{'Most','Least'},...
    'RowNames',{'All Fish','F','P','D'});
writetable(Ftab,[epath 'Hist_Fore_',harv,'_ensem_pset_MostLeast.csv'],...
    'Delimiter',',','WriteRowNames',true)

%Variability and difference
tstats(1,1) = var(tA(:,480));
tstats(2,1) = var(tA(:,1680));
tstats(3,1) = var(tA(:,2880));

tstats(1,2) = max(tA(:,480)) - min(tA(:,480));
tstats(2,2) = max(tA(:,1680)) - min(tA(:,1680));
tstats(3,2) = max(tA(:,2880)) - min(tA(:,2880));

tstats(1,3) = var(tF(:,480));
tstats(2,3) = var(tF(:,1680));
tstats(3,3) = var(tF(:,2880));

tstats(1,4) = max(tF(:,480)) - min(tF(:,480));
tstats(2,4) = max(tF(:,1680)) - min(tF(:,1680));
tstats(3,4) = max(tF(:,2880)) - min(tF(:,2880));

tstats(1,5) = var(tP(:,480));
tstats(2,5) = var(tP(:,1680));
tstats(3,5) = var(tP(:,2880));

tstats(1,6) = max(tP(:,480)) - min(tP(:,480));
tstats(2,6) = max(tP(:,1680)) - min(tP(:,1680));
tstats(3,6) = max(tP(:,2880)) - min(tP(:,2880));

tstats(1,7) = var(tD(:,480));
tstats(2,7) = var(tD(:,1680));
tstats(3,7) = var(tD(:,2880));

tstats(1,8) = max(tD(:,480)) - min(tD(:,480));
tstats(2,8) = max(tD(:,1680)) - min(tD(:,1680));
tstats(3,8) = max(tD(:,2880)) - min(tD(:,2880));

Stab = array2table(tstats,'VariableNames',{'varAll','diffAll','varF','diffF',...
    'varP','diffP','varD','diffD'},...
    'RowNames',{'1900','2000','2100'});
writetable(Stab,[epath 'Hist_Fore_',harv,'_ensem_pset_VarDiff.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% Line color order
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    %0.5 0.5 0.5;...    %med grey
    %49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    %1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

%% moving means
mmtF = movmean(tF,61,2);
mmtP = movmean(tP,61,2);
mmtD = movmean(tD,61,2);
mmtA = movmean(tA,61,2);

mmoF = movmean(tForig,61);
mmoP = movmean(tPorig,61);
mmoD = movmean(tDorig,61);
mmoA = movmean(tAorig,61);

%%
figure(1)
plot(y,log10(mmtA)); hold on;
plot(y,log10(mmoA),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(1) y(end)])
ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','orig')
legend('location','eastoutside')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem.png'])

%%
figure(2)
plot(y,log10(mmtF)); hold on;
plot(y,log10(mmoF),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(1) y(end)])
%ylim([0.425 0.675])
title('Forage fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','orig')
legend('location','eastoutside')
print('-dpng',[ppath 'Hist_Fore_',harv,'_forage_ensem.png'])

%%
figure(3)
plot(y,log10(mmtP)); hold on;
plot(y,log10(mmoP),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(1) y(end)])
%ylim([0.425 0.675])
title('Large pelagic fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','orig')
legend('location','eastoutside')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pel_ensem.png'])

%%
figure(4)
plot(y,log10(mmtD)); hold on;
plot(y,log10(mmoD),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(1) y(end)])
%ylim([0.425 0.675])
title('Demersal fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','orig')
legend('location','eastoutside')
print('-dpng',[ppath 'Hist_Fore_',harv,'_dem_ensem.png'])



