% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3

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
% fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
%     'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];
fpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/';
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat']);

%% Ensemble parameter sets
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);

lpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];

%% ts
%In original saved file
% y1 = 1860+(1/12):(1/12):2005;
% y2 = 2005+(1/12):(1/12):2100;
% y = [y1 y2];

% SOMETHING WRONG WITH MOVING MEAN OF #28 AROUND 1949
% FIX TS OF SMALL FISH AT T=1080
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

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];

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
writetable(Ftab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MostLeast.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ftab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MostLeast.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% Variability and difference
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
writetable(Stab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarDiff.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Stab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarDiff.csv'],...
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
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

%% moving means
mmtF = movmean(tF,31,2); %[1081,28]=1.94e+35
mmtP = movmean(tP,31,2); %[1081,28]=1.94e+35
mmtD = movmean(tD,31,2); %[1081,28]=1.94e+35
mmtA = movmean(tA,31,2); %[1081,28]=5.83e+35

mmoF = movmean(tForig,31);
mmoP = movmean(tPorig,31);
mmoD = movmean(tDorig,31);
mmoA = movmean(tAorig,31);

%% Individual Plots 
% All
figure(1)
plot(y,log10(mmtA)); hold on;
plot(y,log10(mmoA),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1900 2100])
ylim([0.425 0.725])
title('All fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem_mid6_temp3.png'])

%% F
figure(2)
plot(y,log10(mmtF)); hold on;
plot(y,log10(mmoF),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1900 2100])
ylim([0.025 0.325])
title('Forage fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_forage_ensem_mid6_temp3.png'])

%% P
figure(3)
plot(y,log10(mmtP)); hold on;
plot(y,log10(mmoP),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1900 2100])
ylim([-0.275 0.375])
title('Large pelagic fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pel_ensem_mid6_temp3.png'])

%% D
figure(4)
plot(y,log10(mmtD)); hold on;
plot(y,log10(mmoD),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1900 2100])
ylim([-0.15 0.1])
title('Demersal fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_dem_ensem_mid6_temp3.png'])

%% All together on one
figure(5)
% plot(y,log10(mmtA),'k'); hold on;
% plot(y,log10(mmoA),'k','LineWidth',2); hold on;

plot(y,log10(mmtF),'color',[1 80/255 0]); hold on;
plot(y,log10(mmoF),'color',[1 80/255 0],'LineWidth',3); hold on;

plot(y,log10(mmtP),'color',[0 0.5 0.75]); hold on;
plot(y,log10(mmoP),'color',[0 0.5 0.75],'LineWidth',3); hold on;

plot(y,log10(mmtD),'color',[0 0.7 0]); hold on;
plot(y,log10(mmoD),'color',[0 0.5 0],'LineWidth',3); hold on;

xlim([1900 2100])
ylim([-0.3 0.4])
%title('Demersal fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_together_ensem_mid6_temp3.png'])

%% CONE OF UNCERTAINTY
mF = mean(mmtF);
mP = mean(mmtP);
mD = mean(mmtD);
mA = mean(mmtA);

sF = std(mmtF);
sP = std(mmtP);
sD = std(mmtD);
sA = std(mmtA);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%replace with sims with min and max values
Ya=[mmtA(25,:) fliplr(mmtA(12,:))]; 
Yf=[mmtF(24,:) fliplr(mmtF(17,:))]; 
Yp=[mmtP(6,:)  fliplr(mmtP(3,:))]; 
Yd=[mmtD(10,:) fliplr(mmtD(17,:))]; 

%+/- 1 stdev
Sa=[mA+sA fliplr(mA-sA)]; 
Sf=[mF+sF fliplr(mF-sF)]; 
Sp=[mP+sP fliplr(mP-sP)]; 
Sd=[mD+sD fliplr(mD-sD)]; 

%%
figure(6)
f=fill(X,log10(Ya),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,log10(mA),'k','LineWidth',2);
xlim([1900 2100])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem_mid6_temp3_cone_minmax.png'])

%%
figure(7)
fill(X,log10(Yf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,log10(Yp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,log10(Yd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,log10(mF),'color',[1 80/255 0],'LineWidth',2); hold on;
plot(y,log10(mP),'color',[0 0.35 0.75],'LineWidth',2); hold on;
plot(y,log10(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([1900 2100])
ylim([-0.3 0.35])
title('All functional types')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_types_ensem_mid6_temp3_cone_minmax.png'])

%%
figure(8)
f=fill(X,log10(Sa),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,log10(mA),'k','LineWidth',2);
xlim([1900 2100])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem_mid6_temp3_cone_1std.png'])

%%
figure(9)
fill(X,log10(Sf),'r','FaceAlpha',0.25,'EdgeAlpha',0.35); hold on; %plot filled area
fill(X,log10(Sp),'c','FaceAlpha',0.25,'EdgeAlpha',0.35); hold on; %plot filled area
fill(X,log10(Sd),'g','FaceAlpha',0.25,'EdgeAlpha',0.35); hold on; %plot filled area
plot(y,log10(mF),'color',[1 80/255 0],'LineWidth',2); hold on;
plot(y,log10(mP),'color',[0 0.35 0.75],'LineWidth',2); hold on;
plot(y,log10(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([1900 2100])
ylim([-0.2 0.3])
title('All functional types')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_types_ensem_mid6_temp3_cone_1std.png'])



