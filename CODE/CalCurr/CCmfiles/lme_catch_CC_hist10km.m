% Catch on LME scale

clear all
close all

% SAUP
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);     %1950-2006 SAUP average
load([spath 'SAUP_CC_LME_Catch_1988_2010.mat']);    %1988-2010 SAUP avg

%% plot info
%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% 10km grid
hpath = '/Volumes/GFDL/NEMURO/10km/';
load([hpath 'gridspec_10km.mat'],'LON','LAT');
load([hpath 'Data_grid_10km_hist.mat']);

geolon_t = LON;
geolat_t = LAT;
GRD10 = GRD;

% plotminlat=30; %Set these bounds for your data
% plotmaxlat=48;
% plotminlon=-134;
% plotmaxlon=-115;

clear LON LAT GRD

tot_area10_km2 = nansum(GRD10.AREA) * 1e-6;

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
ppath = [pp cfile '/CC/'];

harv = 'All_fish03';
tharv = 'F = 0.3 yr^-^1';

load([fpath 'Means_Historic10km_' harv '_' cfile '.mat'],...
    'mf_my','mp_my','md_my',...
    'lp_my','ld_my');
mf_my10 = mf_my;
mp_my10 = mp_my;
md_my10 = md_my;
lp_my10 = lp_my;
ld_my10 = ld_my;

%% FEISTY LME biomass in MT
% g/m2/d --> total g
Amf_catch10 = mf_my10 .* GRD10.AREA * 365; %mean fish catch per yr
Amp_catch10 = mp_my10 .* GRD10.AREA * 365;
Amd_catch10 = md_my10 .* GRD10.AREA * 365;
Alp_catch10 = lp_my10 .* GRD10.AREA * 365;
Ald_catch10 = ld_my10 .* GRD10.AREA * 365;

%%
AF_catch10 = Amf_catch10;
AP_catch10 = Amp_catch10 + Alp_catch10;
AD_catch10 = Amd_catch10 + Ald_catch10;
All_catch10 = AF_catch10 + AP_catch10 + AD_catch10;

%% g --> MT
AF_mcatch10 = nansum(AF_catch10) * 1e-6;
AP_mcatch10 = nansum(AP_catch10) * 1e-6;
AD_mcatch10 = nansum(AD_catch10) * 1e-6;
All_mcatch10 = nansum(All_catch10) * 1e-6;

%% MT/km2
AF_mcatch10_km2 = AF_mcatch10 ./ tot_area10_km2;
AP_mcatch10_km2 = AP_mcatch10 ./ tot_area10_km2;
AD_mcatch10_km2 = AD_mcatch10 ./ tot_area10_km2;
All_mcatch10_km2 = All_mcatch10 ./ tot_area10_km2;

FracPD10 = AP_mcatch10_km2 ./ (AP_mcatch10_km2 + AD_mcatch10_km2);

%% log-10 transform
l10p_10=log10(All_mcatch10_km2);
l10pF_10=log10(AF_mcatch10_km2);
l10pP_10=log10(AP_mcatch10_km2);
l10pD_10=log10(AD_mcatch10_km2);

%Mean of all years
l10c=log10(nanmean(CC_catch_km2));
l10cF=log10(nanmean(FCC_catch_km2));
l10cP=log10(nanmean(PCC_catch_km2));
l10cD=log10(nanmean(DCC_catch_km2));

%Just CC LME
l10s=log10(slme_mcatch10(3));
l10sF=log10(Flme_mcatch10(3));
l10sP=log10(Plme_mcatch10(3));
l10sD=log10(Dlme_mcatch10(3));

%% Plots
figure(1)
plot(l10sF,l10pF_10,'^k','MarkerSize',10,'MarkerFaceColor','k'); hold on;
plot(l10cF,l10pF_10,'^b','MarkerSize',10,'MarkerFaceColor','b'); hold on;
plot(l10sP,l10pP_10,'vk','MarkerSize',10,'MarkerFaceColor','k'); hold on;
plot(l10cP,l10pP_10,'vb','MarkerSize',10,'MarkerFaceColor','b'); hold on;
plot(l10sD,l10pD_10,'sk','MarkerSize',10,'MarkerFaceColor','k'); hold on;
plot(l10cD,l10pD_10,'sb','MarkerSize',10,'MarkerFaceColor','b'); hold on;
plot(l10s,l10p_10,'ok','MarkerSize',10,'MarkerFaceColor','k'); hold on;
plot(l10c,l10p_10,'ob','MarkerSize',10,'MarkerFaceColor','b'); hold on;
plot(x,x,'--k'); hold on;
plot(x,x2h,':k'); hold on;
plot(x,x2l,':k'); hold on;
plot(x,x5h,':','color',[0.5 0.5 0.5]); hold on;
plot(x,x5l,':','color',[0.5 0.5 0.5]); hold on;
legend('F55','F23','P55','P23','D55','D23','All55','All23')
legend('location','northwest')
axis([-3 1 -3 1])
xlabel('SAU')
ylabel('FEISTY ')
title(['Mean catch (log_1_0 MT km^-^2) with ' tharv])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist10km_',harv,'_SAUP_comp_types.png'])


