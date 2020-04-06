% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Fishing yield in LMEs only 
% total over all months in each year

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
fpath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%% Historic
% load([fpath 'Means_Historic_' harv '_' cfile '.mat'],'time',...
%     'mf_tsyc','mp_tsyc','md_tsyc','lp_tsyc','ld_tsyc',...
%     'units_yield','units_catch');
% 
% HF_tsyc = mf_tsyc;
% HP_tsyc = mp_tsyc + lp_tsyc;
% HD_tsyc = md_tsyc + ld_tsyc;
% HA_tsyc = HF_tsyc + HP_tsyc + HD_tsyc;
% 
% HMF_tsyc = mf_tsyc;
% HLP_tsyc = lp_tsyc;
% HLD_tsyc = ld_tsyc;
% HAA_tsyc = HMF_tsyc + HLP_tsyc + HLD_tsyc;
% 
% clear mf_tsyc mp_tsyc md_tsyc lp_tsyc ld_tsyc

%% Forecast
% load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'time',...
%     'mf_tsyc','mp_tsyc','md_tsyc','lp_tsyc','ld_tsyc',...
%     'units_yield','units_catch');
% 
% FF_tsyc = mf_tsyc;
% FP_tsyc = mp_tsyc + lp_tsyc;
% FD_tsyc = md_tsyc + ld_tsyc;
% FA_tsyc = FF_tsyc + FP_tsyc + FD_tsyc;
% 
% FMF_tsyc = mf_tsyc;
% FLP_tsyc = lp_tsyc;
% FLD_tsyc = ld_tsyc;
% FAA_tsyc = FMF_tsyc + FLP_tsyc + FLD_tsyc;
% 
% clear mf_tsyc mp_tsyc md_tsyc lp_tsyc ld_tsyc
% 
% save([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
%     'HF_tsyc','HP_tsyc','HD_tsyc','HA_tsyc','HMF_tsyc','HLP_tsyc',...
%     'HLD_tsyc','HAA_tsyc',...
%     'FF_tsyc','FP_tsyc','FD_tsyc','FA_tsyc','FMF_tsyc','FLP_tsyc',...
%     'FLD_tsyc','FAA_tsyc','-append');
% 
% save([fpath2 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
%     'HF_tsyc','HP_tsyc','HD_tsyc','HA_tsyc','HMF_tsyc','HLP_tsyc',...
%     'HLD_tsyc','HAA_tsyc',...
%     'FF_tsyc','FP_tsyc','FD_tsyc','FA_tsyc','FMF_tsyc','FLP_tsyc',...
%     'FLD_tsyc','FAA_tsyc','-append');

load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
    'HF_tsyc','HP_tsyc','HD_tsyc','HA_tsyc','HMF_tsyc','HLP_tsyc',...
    'HLD_tsyc','HAA_tsyc','FF_tsyc','FP_tsyc','FD_tsyc','FA_tsyc',...
    'FMF_tsyc','FLP_tsyc','FLD_tsyc','FAA_tsyc');

%% Ensemble parameter sets
efile = 'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';

epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
epath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];

load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Hlme_tsc_mf','Hlme_tsc_mp','Hlme_tsc_md','Hlme_tsc_lp','Hlme_tsc_ld');
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Flme_tsc_mf','Flme_tsc_mp','Flme_tsc_md','Flme_tsc_lp','Flme_tsc_ld');

% Time series of catch (g per km2 per year)

%% ts
%In original saved file
y1 = 1861:2005;
y2 = 2006:2100;
y = [y1 y2];

HF = Hlme_tsc_mf;
HP = Hlme_tsc_mp + Hlme_tsc_lp;
HD = Hlme_tsc_md + Hlme_tsc_ld;
HLP = Hlme_tsc_lp;
HLD = Hlme_tsc_ld;

%% SOMETHING WRONG WITH MEAN OF #28 AROUND 1950
% FIX TS AT T=90 & 91 w/interp1
yi = y(89:92);
HF(28,89:92) = interp1(y([89,92]),HF(28,[89,92]),yi);
HP(28,89:92) = interp1(y([89,92]),HP(28,[89,92]),yi);
HD(28,89:92) = interp1(y([89,92]),HD(28,[89,92]),yi);
HLP(28,89:92) = interp1(y([89,92]),HLP(28,[89,92]),yi);
HLD(28,89:92) = interp1(y([89,92]),HLD(28,[89,92]),yi);

%%
HA = HF + HP + HD;
HAA = HF + HLP + HLD;

FF = Flme_tsc_mf;
FP = Flme_tsc_mp + Flme_tsc_lp;
FD = Flme_tsc_md + Flme_tsc_ld;
FA = FF + FP + FD;
FLP = Flme_tsc_lp;
FLD = Flme_tsc_ld;
FAA = FF + FLP + FLD;

tForig = [HF_tsyc FF_tsyc];
tPorig = [HP_tsyc FP_tsyc];
tDorig = [HD_tsyc FD_tsyc];
tAorig = [HA_tsyc FA_tsyc];
tLPorig = [HLP_tsyc FLP_tsyc];
tLDorig = [HLD_tsyc FLD_tsyc];
tAAorig = [HAA_tsyc FAA_tsyc];

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];
tLP = [HLP FLP];
tLD = [HLD FLD];
tAA = [HAA FAA];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];
tLP = [tLP; tLPorig];
tLD = [tLD; tLDorig];
tAA = [tAA; tAAorig];

%% difference from 1951
test=find(y>1950);
yid=test(1);

dtF = tF - tF(:,yid);
dtP = tP - tP(:,yid);
dtD = tD - tD(:,yid);
dtA = tA - tA(:,yid);
dtLP = tLP - tLP(:,yid);
dtLD = tLD - tLD(:,yid);
dtAA = tAA - tAA(:,yid);

% dtF = tForig - tForig(:,yid);
% dtP = tPorig - tPorig(:,yid);
% dtD = tDorig - tDorig(:,yid);
% dtA = tAorig - tAorig(:,yid);

% percent difference from 1951
pdF = (tF - tF(:,yid)) ./ tF(:,yid); 
pdP = (tP - tP(:,yid)) ./ tP(:,yid);
pdD = (tD - tD(:,yid)) ./ tD(:,yid); 
pdA = (tA - tA(:,yid)) ./ tA(:,yid);
pdLP = (tLP - tLP(:,yid)) ./ tLP(:,yid);
pdLD = (tLD - tLD(:,yid)) ./ tLD(:,yid); 
pdAA = (tAA - tAA(:,yid)) ./ tAA(:,yid);

%% Cone of uncert raw diff
rmF = mean(dtF);
rmP = mean(dtP);
rmD = mean(dtD);
rmA = mean(dtA);
rmLP = mean(dtLP);
rmLD = mean(dtLD);
rmAA = mean(dtAA);

rsF = std(dtF);
rsP = std(dtP);
rsD = std(dtD);
rsA = std(dtA);
rsLP = std(dtLP);
rsLD = std(dtLD);
rsAA = std(dtAA);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%+/- 1 stdev
Ra=[rmA+rsA fliplr(rmA-rsA)]; 
Rf=[rmF+rsF fliplr(rmF-rsF)]; 
Rp=[rmP+rsP fliplr(rmP-rsP)]; 
Rd=[rmD+rsD fliplr(rmD-rsD)];
Raa=[rmAA+rsAA fliplr(rmAA-rsAA)]; 
Rlp=[rmLP+rsLP fliplr(rmLP-rsLP)]; 
Rld=[rmLD+rsLD fliplr(rmLD-rsLD)]; 

%% Cone pdiff
mpF = mean(pdF);
mpP = mean(pdP);
mpD = mean(pdD);
mpA = mean(pdA);
mpLP = mean(pdLP);
mpLD = mean(pdLD);
mpAA = mean(pdAA);

spF = std(pdF);
spP = std(pdP);
spD = std(pdD);
spA = std(pdA);
spLP = std(pdLP);
spLD = std(pdLD);
spAA = std(pdAA);

%create y values for out and then back
%+/- 1 stdev
Va=[mpA+spA fliplr(mpA-spA)]; 
Vf=[mpF+spF fliplr(mpF-spF)]; 
Vp=[mpP+spP fliplr(mpP-spP)]; 
Vd=[mpD+spD fliplr(mpD-spD)]; 
Vaa=[mpAA+spAA fliplr(mpAA-spAA)]; 
Vlp=[mpLP+spLP fliplr(mpLP-spLP)]; 
Vld=[mpLD+spLD fliplr(mpLD-spLD)]; 

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
%Difference 
(tForig(end)-tForig(:,yid))/tForig(:,yid)
(tDorig(end)-tDorig(:,yid))/tDorig(:,yid)
(tPorig(end)-tPorig(:,yid))/tPorig(:,yid)

%% all - baseline diff
figure(1)
subplot(2,2,1)
plot(y,(tA(end,:)*1e-6),'k','LineWidth',2); hold on;
xlim([y(yid) y(end)])
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1)')
title('All fishes')

% types - baseline diff
subplot(2,2,2)
plot(y,(tF(end,:)*1e-6),'r','LineWidth',2); hold on;
plot(y,(tP(end,:)*1e-6),'b','LineWidth',2); hold on;
plot(y,(tD(end,:)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')

subplot(2,2,3)
plot(y,(dtA(end,:)*1e-6),'k','LineWidth',2); hold on;
xlim([y(yid) y(end)])
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')

% types - baseline diff
subplot(2,2,4)
plot(y,(dtF(end,:)*1e-6),'r','LineWidth',2); hold on;
plot(y,(dtP(end,:)*1e-6),'b','LineWidth',2); hold on;
plot(y,(dtD(end,:)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_types_subplot.png'])

%% baseline adults diff
figure(2)
subplot(2,2,1)
plot(y,(tAA(end,:)*1e-6),'k','LineWidth',2); hold on;
xlim([y(yid) y(end)])
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1)')
title('Adult fishes')

% types - baseline diff
subplot(2,2,2)
plot(y,(tF(end,:)*1e-6),'r','LineWidth',2); hold on;
plot(y,(tLP(end,:)*1e-6),'b','LineWidth',2); hold on;
plot(y,(tLD(end,:)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')

subplot(2,2,3)
plot(y,(dtAA(end,:)*1e-6),'k','LineWidth',2); hold on;
xlim([y(yid) y(end)])
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')

% types - baseline diff
subplot(2,2,4)
plot(y,(dtF(end,:)*1e-6),'r','LineWidth',2); hold on;
plot(y,(dtLP(end,:)*1e-6),'b','LineWidth',2); hold on;
plot(y,(dtLD(end,:)*1e-6),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_adults_subplot.png'])

%% types - baseline pdiff
figure(3)
subplot(2,2,1)
plot(y,(pdA(end,:)),'k','LineWidth',2); hold on;
xlim([y(yid) y(end)])
xlabel('Year')
ylabel('% \Delta in yield (MT km^-^2 y^-^1) relative to 1951')
title('All fishes')

% types - baseline diff
subplot(2,2,2)
plot(y,(pdF(end,:)),'r','LineWidth',2); hold on;
plot(y,(pdP(end,:)),'b','LineWidth',2); hold on;
plot(y,(pdD(end,:)),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_yield_types_subplot.png'])

%% baseline adults pdiff
figure(4)
subplot(2,2,1)
plot(y,(pdAA(end,:)),'k','LineWidth',2); hold on;
xlim([y(yid) y(end)])
xlabel('Year')
ylabel('% \Delta in yield (MT km^-^2 y^-^1) relative to 1951')
title('Adult fishes')

% types - baseline diff
subplot(2,2,2)
plot(y,(pdF(end,:)),'r','LineWidth',2); hold on;
plot(y,(pdLP(end,:)),'b','LineWidth',2); hold on;
plot(y,(pdLD(end,:)),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
%title('All functional types')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_yield_adults_subplot.png'])

%% all w/uncert
figure(5)
f=fill(X,(Ra)*1e-6,'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(rmA)*1e-6,'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('All fishes')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_all_ensem_mid6_temp3_cone_1std_yr.png'])

%% adults w/uncert
figure(6)
f=fill(X,(Raa)*1e-6,'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(rmAA)*1e-6,'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('Adult fishes')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_adult_ensem_mid6_temp3_cone_1std_yr.png'])

%% all w/uncert - pdiff
figure(7)
f=fill(X,(Va),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(mpA),'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('All fish')
xlabel('Year')
ylabel('Percent change in yield relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_yield_all_ensem_mid6_temp3_cone_1std_yr.png'])

%% adults w/uncert - pdiff
figure(8)
f=fill(X,(Vaa),'k','FaceAlpha',0.25,'EdgeAlpha',0.25);  %plot filled area
hold on
plot(y,(mpAA),'k','LineWidth',2);
xlim([y(yid) y(end)])
% ylim([0.425 0.675])
title('Adult fish')
xlabel('Year')
ylabel('Percent change in yield relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_yield_adult_ensem_mid6_temp3_cone_1std_yr.png'])

%% types - diff
figure(9)
fill(X,(Rf)*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rp)*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rd)*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(rmF)*1e-6,'r','LineWidth',2); hold on;
plot(y,(rmP)*1e-6,'b','LineWidth',2); hold on;
plot(y,(rmD)*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
title('All functional types and stages')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_types_all_ensem_mid6_temp3_cone_1std_yr_v2.png'])

%% adults - diff
figure(10)
fill(X,(Rf)*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rlp)*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rld)*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(rmF)*1e-6,'r','LineWidth',2); hold on;
plot(y,(rmLP)*1e-6,'b','LineWidth',2); hold on;
plot(y,(rmLD)*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
title('Adults of all functional types')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_diff1951_yield_types_adults_ensem_mid6_temp3_cone_1std_yr_v2.png'])

%% types - pdiff
figure(11)
fill(X,(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(mpF),'r','LineWidth',2); hold on;
plot(y,(mpP),'b','LineWidth',2); hold on;
plot(y,(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
title('All functional types and stages')
xlabel('Year')
ylabel('Percent change in yield relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_yield_types_all_ensem_mid6_temp3_cone_1std_yr.png'])

%% adults - pdiff
figure(12)
fill(X,(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,(mpF),'r','LineWidth',2); hold on;
plot(y,(mpP),'b','LineWidth',2); hold on;
plot(y,(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
title('Adults of all functional types')
xlabel('Year')
ylabel('Percent change in yield relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pdiff1951_yield_types_adults_ensem_mid6_temp3_cone_1std_yr.png'])

%% save 
save([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_yield.mat'],...
    'Rf','Rp','Rd','Ra','Rlp','Rld','Raa',...
    'rmF','rmP','rmD','rmA','rmLP','rmLD','rmAA',...
    'dtF','dtP','dtD','dtA','dtLP','dtLD','dtAA',...
    'Vf','Vp','Vd','Va','Vlp','Vld','Vaa',...
    'mpF','mpP','mpD','mpA','mpLP','mpLD','mpAA','X','y');
save([epath2 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_yield.mat'],...
    'Rf','Rp','Rd','Ra','Rlp','Rld','Raa',...
    'rmF','rmP','rmD','rmA','rmLP','rmLD','rmAA',...
    'dtF','dtP','dtD','dtA','dtLP','dtLD','dtAA',...
    'Vf','Vp','Vd','Va','Vlp','Vld','Vaa',...
    'mpF','mpP','mpD','mpA','mpLP','mpLD','mpAA','X','y');



