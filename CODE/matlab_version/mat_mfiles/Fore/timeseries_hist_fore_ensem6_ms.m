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
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat']);

%% Ensemble parameter sets
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);

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
% set(groot,'defaultFontName','TimesNewRoman');
% set(groot,'defaultFontSize',16);

%% moving means
mmtF = movmean(tF,31,2); %[1081,28]=1.94e+35
mmtP = movmean(tP,31,2); %[1081,28]=1.94e+35
mmtD = movmean(tD,31,2); %[1081,28]=1.94e+35
mmtA = movmean(tA,31,2); %[1081,28]=5.83e+35

mmoF = movmean(tForig,31);
mmoP = movmean(tPorig,31);
mmoD = movmean(tDorig,31);
mmoA = movmean(tAorig,31);

%% SubPlots of ensemble of all types
% All
figure(1)
subplot(2,2,4)
plot(y,log10(mmtA)); hold on;
plot(y,log10(mmoA),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1951 2100])
ylim([0.425 0.725])
title('All fish')
xlabel('Year')

% F
subplot(2,2,1)
plot(y,log10(mmtF)); hold on;
plot(y,log10(mmoF),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1951 2100])
ylim([0.025 0.325])
title('Forage fish')
ylabel('log10 Biomass (g m^-^2)')

% P
subplot(2,2,2)
plot(y,log10(mmtP)); hold on;
plot(y,log10(mmoP),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1951 2100])
ylim([-0.275 0.375])
title('Large pelagic fish')

% D
subplot(2,2,3)
plot(y,log10(mmtD)); hold on;
plot(y,log10(mmoD),'color',[0 0.5 0.75],'LineWidth',2);
xlim([1951 2100])
ylim([-0.15 0.1])
title('Demersal fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_subplot_all_ensem_mid6_temp3.png'])

%% Orig params All together on one
figure(2)
plot(y,log10(tForig),'r','LineWidth',2); hold on;
plot(y,log10(tDorig),'color',[0 0.7 0],'LineWidth',2); hold on;
plot(y,log10(tPorig),'b','LineWidth',2); hold on;
legend('Forage','Demersal','Large Pelagic')
xlim([1951 2100])
ylim([-0.15 0.3])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_together_orig_param_legend.png'])

figure(3)
plot(y,log10(tForig),'r','LineWidth',2); hold on;
plot(y,log10(tDorig),'color',[0 0.7 0],'LineWidth',2); hold on;
plot(y,log10(tPorig),'b','LineWidth',2); hold on;
xlim([1951 2100])
ylim([-0.15 0.3])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_together_orig_param_noleg.png'])

%% CONE OF UNCERTAINTY
mF = mean(tF);
mP = mean(tP);
mD = mean(tD);
mA = mean(tA);

sF = std(mmtF);
sP = std(mmtP);
sD = std(mmtD);
sA = std(mmtA);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%replace with sims with min and max values
Ya=[tA(25,:) fliplr(tA(12,:))]; 
Yf=[tF(24,:) fliplr(tF(17,:))]; 
Yp=[tP(6,:)  fliplr(tP(3,:))]; 
Yd=[tD(10,:) fliplr(tD(17,:))]; 

%+/- 1 stdev
Sa=[mA+sA fliplr(mA-sA)]; 
Sf=[mF+sF fliplr(mF-sF)]; 
Sp=[mP+sP fliplr(mP-sP)]; 
Sd=[mD+sD fliplr(mD-sD)]; 

%%
figure(4)
f=fill(X,log10(Ya),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,log10(mA),'k','LineWidth',2);
xlim([1951 2100])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem_mid6_temp3_cone_minmax_ms.png'])

%%
figure(5)
fill(X,log10(Yf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,log10(Yp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,log10(Yd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,log10(mF),'r','LineWidth',2); hold on;
plot(y,log10(mP),'b','LineWidth',2); hold on;
plot(y,log10(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([1951 2100])
ylim([-0.3 0.35])
%title('All functional types')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_types_ensem_mid6_temp3_cone_minmax_ms.png'])

%%
figure(6)
f=fill(X,log10(Sa),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,log10(mA),'k','LineWidth',2);
xlim([1951 2100])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem_mid6_temp3_cone_1std_ms.png'])

%%
figure(7)
fill(X,log10(Sf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,log10(Sp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,log10(Sd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,log10(mF),'r','LineWidth',2); hold on;
plot(y,log10(mP),'b','LineWidth',2); hold on;
plot(y,log10(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([1951 2100])
ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_types_ensem_mid6_temp3_cone_1std_ms.png'])



