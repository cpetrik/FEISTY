% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Production instead of biomass

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
%epath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%% Historic
% load([fpath 'Means_Historic_',harv,'_prod_' cfile '.mat'],...
%     'sf_tamean','sp_tamean','sd_tamean',...
%     'mf_tamean','mp_tamean','md_tamean',...
%     'lp_tamean','ld_tamean','b_tamean');
% 
% HF_tamean = sf_tamean + mf_tamean;
% HP_tamean = sp_tamean + mp_tamean + lp_tamean;
% HD_tamean = sd_tamean + md_tamean + ld_tamean;
% HB_tamean = b_tamean;
% HA_tamean = HF_tamean + HP_tamean + HD_tamean;
% 
% clear sf_tamean sp_tamean sd_tamean mf_tamean mp_tamean md_tamean lp_tamean ld_tamean b_tamean

%% Forecast
% load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
%     'sf_tamean','sp_tamean','sd_tamean',...
%     'mf_tamean','mp_tamean','md_tamean',...
%     'lp_tamean','ld_tamean','b_tamean');
% 
% FF_tamean = sf_tamean + mf_tamean;
% FP_tamean = sp_tamean + mp_tamean + lp_tamean;
% FD_tamean = sd_tamean + md_tamean + ld_tamean;
% FB_tamean = b_tamean;
% FA_tamean = FF_tamean + FP_tamean + FD_tamean;
% 
% clear sf_tamean sp_tamean sd_tamean mf_tamean mp_tamean md_tamean lp_tamean ld_tamean b_tamean
% 
% save([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
%     'HF_tamean','HP_tamean','HD_tamean','HB_tamean',...
%     'FF_tamean','FP_tamean','FD_tamean','FB_tamean',...
%     'HA_tamean','FA_tamean','-append');
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
    'HF_tamean','HP_tamean','HD_tamean','HB_tamean',...
    'FF_tamean','FP_tamean','FD_tamean','FB_tamean',...
    'HA_tamean','FA_tamean');

%% Ensemble parameter sets
% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'haTsF','haTsP','haTsD','haTmF','haTmP','haTmD','haTB','haTlP','haTlD');
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'faTsF','faTsP','faTsD','faTmF','faTmP','faTmD','faTB','faTlP','faTlD');

%aT=tamean=nansum(prod.*area,1);

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

%% Prod in g/d --> g/yr?
tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];

%% difference from 1951
test=find(y>1950);
yid=test(1);

dtF = tF - tF(:,yid);
dtP = tP - tP(:,yid);
dtD = tD - tD(:,yid);
dtA = tA - tA(:,yid);

% dtF = tForig - tForig(:,yid);
% dtP = tPorig - tPorig(:,yid);
% dtD = tDorig - tDorig(:,yid);
% dtA = tAorig - tAorig(:,yid);

% percent difference from 1951
pdF = (tF - tF(:,yid)) ./ tF(:,yid); 
pdP = (tP - tP(:,yid)) ./ tP(:,yid);
pdD = (tD - tD(:,yid)) ./ tD(:,yid); 
pdA = (tA - tA(:,yid)) ./ tA(:,yid);

%% moving means
mmtF = movmean(tF,12,2); 
mmtP = movmean(tP,12,2); 
mmtD = movmean(tD,12,2); 
mmtA = movmean(tA,12,2); 

mmoF = movmean(tForig,12);
mmoP = movmean(tPorig,12);
mmoD = movmean(tDorig,12);
mmoA = movmean(tAorig,12);

% difference from 1951
mmdF = mmtF - mmtF(:,yid); 
mmdP = mmtP - mmtP(:,yid);
mmdD = mmtD - mmtD(:,yid); 
mmdA = mmtA - mmtA(:,yid);

moF = mmoF - mmoF(:,yid); 
moP = mmoP - mmoP(:,yid);
moD = mmoD - mmoD(:,yid); 
moA = mmoA - mmoA(:,yid);

% percent difference from 1951
mpdF = (mmtF - mmtF(:,yid)) ./ mmtF(:,yid); 
mpdP = (mmtP - mmtP(:,yid)) ./ mmtP(:,yid);
mpdD = (mmtD - mmtD(:,yid)) ./ mmtD(:,yid); 
mpdA = (mmtA - mmtA(:,yid)) ./ mmtA(:,yid);
% mpdF = movmean(pdF,12,2); 
% mpdP = movmean(pdP,12,2); 
% mpdD = movmean(pdD,12,2); 
% mpdA = movmean(pdA,12,2); 

%% Cone of uncert raw diff
rmF = mean(dtF);
rmP = mean(dtP);
rmD = mean(dtD);
rmA = mean(dtA);

rsF = std(dtF);
rsP = std(dtP);
rsD = std(dtD);
rsA = std(dtA);

%create y values for out and then back
%+/- 1 stdev
Ra=[rmA+rsA fliplr(rmA-rsA)]; 
Rf=[rmF+rsF fliplr(rmF-rsF)]; 
Rp=[rmP+rsP fliplr(rmP-rsP)]; 
Rd=[rmD+rsD fliplr(rmD-rsD)]; 

%% CONE OF UNCERTAINTY moving means
mF = mean(mmdF);
mP = mean(mmdP);
mD = mean(mmdD);
mA = mean(mmdA);

% sF = std(mmdF);
% sP = std(mmdP);
% sD = std(mmdD);
% sA = std(mmdA);
sF = movmean(rsF,12,2); 
sP = movmean(rsP,12,2); 
sD = movmean(rsD,12,2); 
sA = movmean(rsA,12,2); 

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%+/- 1 stdev
Sa=[mA+sA fliplr(mA-sA)]; 
Sf=[mF+sF fliplr(mF-sF)]; 
Sp=[mP+sP fliplr(mP-sP)]; 
Sd=[mD+sD fliplr(mD-sD)]; 

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

%% 
% NOT THE SAME AS TROPHIC AMP BAR GRAPH
% MAY NEED TOT(PROD.*AREA) INSTEAD OF MEAN(PROD.*AREA)
% Difference from Troph Amp seems to be related to pure diff vs. % diff
(tForig(end)-tForig(:,yid))/tForig(:,yid)
(tDorig(end)-tDorig(:,yid))/tDorig(:,yid)
(tPorig(end)-tPorig(:,yid))/tPorig(:,yid)

%% types - baseline diff
figure(1)
plot(y,(moF),'r','LineWidth',2); hold on;
plot(y,(moP),'b','LineWidth',2); hold on;
plot(y,(moD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel('Production (g d^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_prod_types.png'])

%% types - baseline pdiff
figure(2)
plot(y,(mpdF(end,:)),'r','LineWidth',2); hold on;
plot(y,(mpdP(end,:)),'b','LineWidth',2); hold on;
plot(y,(mpdD(end,:)),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
ylabel('Percent change in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_prod_types.png'])

%% all w/uncert
figure(3)
f=fill(X,(Sa),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(mA),'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('Production (g d^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_prod_all_ensem_mid6_temp3_cone_1std_yr.png'])

%% all w/uncert - pdiff
figure(4)
f=fill(X,(Va),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(mpA),'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('Percent change in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_prod_all_ensem_mid6_temp3_cone_1std_yr.png'])

%% types - diff
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
ylabel('Production (g d^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_prod_types_ensem_mid6_temp3_cone_1std_yr_v2.png'])

%% types - pdiff
figure(6)
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
ylabel('Percent change in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_prod_types_ensem_mid6_temp3_cone_1std_yr.png'])

%% save for multipanel plot w/biom changes
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
save([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod.mat'],...
    'Sf','Sp','Sd','Sa','mF','mP','mD','mA','dtF','dtP','dtD','dtA',...
    'Vf','Vp','Vd','Va','mpF','mpP','mpD','mpA','X','y');
save([dpath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod.mat'],...
    'Sf','Sp','Sd','Sa','mF','mP','mD','mA','dtF','dtP','dtD','dtA',...
    'Vf','Vp','Vd','Va','mpF','mpP','mpD','mpA','X','y');



