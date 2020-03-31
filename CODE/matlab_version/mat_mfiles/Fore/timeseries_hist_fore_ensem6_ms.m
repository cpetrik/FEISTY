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
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
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

% difference from 1951
test=find(y>1950);
yid=test(1);
dtF = tF - tF(:,yid);
dtP = tP - tP(:,yid);
dtD = tD - tD(:,yid);
dtA = tA - tA(:,yid);

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
mmtF = movmean(tF,12,2); %[1081,28]=1.94e+35
mmtP = movmean(tP,12,2); %[1081,28]=1.94e+35
mmtD = movmean(tD,12,2); %[1081,28]=1.94e+35
mmtA = movmean(tA,12,2); %[1081,28]=5.83e+35

mmoF = movmean(tForig,12);
mmoP = movmean(tPorig,12);
mmoD = movmean(tDorig,12);
mmoA = movmean(tAorig,12);

% difference from 1951
mmdF = mmtF - mmtF(:,yid); %[1081,28]=1.94e+35
mmdP = mmtP - mmtP(:,yid); %[1081,28]=1.94e+35
mmdD = mmtD - mmtD(:,yid); %[1081,28]=1.94e+35
mmdA = mmtA - mmtA(:,yid);

% percent difference from 1951
mpdF = (mmtF - mmtF(:,yid)) ./ mmtF(:,yid); 
mpdP = (mmtP - mmtP(:,yid)) ./ mmtP(:,yid);
mpdD = (mmtD - mmtD(:,yid)) ./ mmtD(:,yid); 
mpdA = (mmtA - mmtA(:,yid)) ./ mmtA(:,yid);

%% SubPlots of ensemble of all types
% All
figure(1)
subplot(2,2,4)
plot(y,(mmdA)); hold on;
%plot(y,(mmdA(end,:)),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(yid) y(end)])
%ylim([0.425 0.725])
title('All fish')
xlabel('Year')

% F
subplot(2,2,1)
plot(y,(mmdF)); hold on;
%plot(y,(mmdF(end,:)),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(yid) y(end)])
%ylim([0.025 0.325])
title('Forage fish')
ylabel(' Biomass (g m^-^2)')

% P
subplot(2,2,2)
plot(y,(mmdP)); hold on;
%plot(y,(mmdP(end,:)),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(yid) y(end)])
%ylim([-0.275 0.375])
title('Large pelagic fish')

% D
subplot(2,2,3)
plot(y,(mmdD)); hold on;
%plot(y,(mmdD(end,:)),'color',[0 0.5 0.75],'LineWidth',2);
xlim([y(yid) y(end)])
%ylim([-0.15 0.1])
title('Demersal fish')
xlabel('Year')
ylabel('Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_subplot_diff1951_all_ensem_mid6_temp3.png'])

%% Orig params All together on one
figure(2)
plot(y,(dtF(end,:)),'r','LineWidth',2); hold on;
plot(y,(dtD(end,:)),'color',[0 0.7 0],'LineWidth',2); hold on;
plot(y,(dtP(end,:)),'b','LineWidth',2); hold on;
legend('Forage','Demersal','Large Pelagic')
xlim([y(yid) y(end)])
%ylim([-0.15 0.3])
xlabel('Year')
ylabel('Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_all_together_orig_param_legend.png'])

figure(3)
plot(y,(dtF(end,:)),'r','LineWidth',2); hold on;
plot(y,(dtD(end,:)),'color',[0 0.7 0],'LineWidth',2); hold on;
plot(y,(dtP(end,:)),'b','LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.15 0.3])
xlabel('Year')
ylabel('Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_all_together_orig_param_noleg.png'])

%% CONE OF UNCERTAINTY Raw
mF = mean(dtF);
mP = mean(dtP);
mD = mean(dtD);
mA = mean(dtA);

sF = std(dtF);
sP = std(dtP);
sD = std(dtD);
sA = std(dtA);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%+/- 1 stdev
Sa=[mA+sA fliplr(mA-sA)]; 
Sf=[mF+sF fliplr(mF-sF)]; 
Sp=[mP+sP fliplr(mP-sP)]; 
Sd=[mD+sD fliplr(mD-sD)]; 

%% all
figure(4)
f=fill(X,(Sa),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(mA),'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel(' Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_all_types_ensem_mid6_temp3_cone_1std_month.png'])

%% types
figure(5)
fill(X,(Sf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Sp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Sd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(mF),'r','LineWidth',2); hold on;
plot(y,(mP),'b','LineWidth',2); hold on;
plot(y,(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel(' Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_types_ensem_mid6_temp3_cone_1std_month.png'])


%% CONE OF UNCERTAINTY moving means
mF = mean(mmdF);
mP = mean(mmdP);
mD = mean(mmdD);
mA = mean(mmdA);

% msF = std(mmdF);
% msP = std(mmdP);
% msD = std(mmdD);
% msA = std(mmdA);
msF = movmean(sF,12,2); 
msP = movmean(sP,12,2); 
msD = movmean(sD,12,2); 
msA = movmean(sA,12,2);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%+/- 1 stdev
Sa=[mA+msA fliplr(mA-msA)]; 
Sf=[mF+msF fliplr(mF-msF)]; 
Sp=[mP+msP fliplr(mP-msP)]; 
Sd=[mD+msD fliplr(mD-msD)]; 

%% all
figure(6)
f=fill(X,(Sa),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(mA),'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('Biomass (g m^-^2) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_all_types_ensem_mid6_temp3_cone_1std_yr.png'])

%% types
figure(7)
fill(X,(Sf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Sp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Sd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(mF),'r','LineWidth',2); hold on;
plot(y,(mP),'b','LineWidth',2); hold on;
plot(y,(mD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel('Biomass (g m^-^2) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_types_ensem_mid6_temp3_cone_1std_yr_v2.png'])

%% Cone moving means pdiff
mpF = mean(mpdF);
mpP = mean(mpdP);
mpD = mean(mpdD);
mpA = mean(mpdA);

spF = std(mpdF);
spP = std(mpdP);
spD = std(mpdD);
spA = std(mpdA);

%create y values for out and then back
%+/- 1 stdev
Va=[mpA+spA fliplr(mpA-spA)]; 
Vf=[mpF+spF fliplr(mpF-spF)]; 
Vp=[mpP+spP fliplr(mpP-spP)]; 
Vd=[mpD+spD fliplr(mpD-spD)]; 

%% types - pdiff
figure(8)
fill(X,(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(mpF),'r','LineWidth',2); hold on;
plot(y,(mpP),'b','LineWidth',2); hold on;
plot(y,(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel('Percent change in biomass relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_biom_types_ensem_mid6_temp3_cone_1std_yr_v2.png'])


%% save for multipanel plot w/prod changes
save([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_biomass.mat'],...
    'Sf','Sp','Sd','Sa','mF','mP','mD','mA','dtF','dtP','dtD','dtA',...
    'Vf','Vp','Vd','Va','mpF','mpP','mpD','mpA');
save([fpath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_biomass.mat'],...
    'Sf','Sp','Sd','Sa','mF','mP','mD','mA','dtF','dtP','dtD','dtA',...
    'Vf','Vp','Vd','Va','mpF','mpP','mpD','mpA');

