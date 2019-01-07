%POEM catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

AREA_OCN = max(area,1);

% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);
lme_area_km2 = lme_area * 1e-6;


%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

cmYOR=cbrewer('seq','YlOrRd',28);
cmRP=cbrewer('seq','RdPu',28);
cmPR=cbrewer('seq','PuRd',28);

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP MT/km2
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
plme_Bmbio = lme_mbio(:,9) * 1e-6;
plme_Bsbio = lme_sbio(:,9) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
plme_Bmbio = plme_Bmbio ./ lme_area_km2;
plme_Bsbio = plme_Bsbio ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
% Stats
%r
rall=corr(l10s(keep),l10p(keep));
rF=corr(l10sF(keep),l10pF(keep));
rP=corr(l10sP(keep),l10pP(keep));
rD=corr(l10sD(keep),l10pD(keep));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

%% Plots
% All LMEs
figure(1)
plot(x,x,'--k'); hold on;
%scatter(l10s,l10p,30,lme_ptemp(:,1),'filled'); hold on;
scatter(l10s,l10p,30,'k','filled'); hold on;
% cmocean('thermal');
% colorbar('location','eastoutside')
% caxis([-2 28])
axis([-2 2 -2 2])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_all_temp_Stock_ppt.png'])


% no LELC
figure(2)
plot(x,x,'--k'); hold on;
%scatter(l10s(keep),l10p(keep),30,lme_ptemp(keep,1),'filled'); hold on;
scatter(l10s(keep),l10p(keep),30,'k','filled'); hold on;
% cmocean('thermal');
text(-1.75,1.4,['r = ' sprintf('%2.2f',rall)])
text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 2 -2 2])
% colorbar('location','eastoutside')
% caxis([-2 28])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_all_temp_Stock_LELC_ppt.png'])


%% For ms
figure(3)
subplot(2,2,2)
plot(x,x,'--k'); hold on;
%scatter(l10sF(keep),l10pF(keep),25,lme_ptemp(keep,1),'filled'); hold on;
scatter(l10sF(keep),l10pF(keep),25,'k','filled'); hold on;
% cmocean('thermal');
% caxis([-2 28])
%colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-6 2 -6 2])
title('Forage Fishes')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
%scatter(l10sP(keep),l10pP(keep),25,lme_ptemp(keep,1),'filled'); hold on;
scatter(l10sP(keep),l10pP(keep),25,'k','filled'); hold on;
% cmocean('thermal');
% caxis([-2 28])
text(-5.5,1.5,['r = ' sprintf('%2.2f',rP)])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 2 -6 2])
title('Large Pelagics')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
%scatter(l10sD(keep),l10pD(keep),25,lme_ptemp(keep,1),'filled'); hold on;
scatter(l10sD(keep),l10pD(keep),25,'k','filled'); hold on;
% cmocean('thermal');
% caxis([-2 28])
text(-1.75,1.7,['r = ' sprintf('%2.2f',rD)])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 2 -2 2])
title('Demersals')

subplot(2,2,1)
plot(x,x,'--k'); hold on;
%scatter(l10s(keep),l10p(keep),25,lme_ptemp(keep,1),'filled'); hold on;
scatter(l10s(keep),l10p(keep),25,'k','filled'); hold on;
% cmocean('thermal');
% caxis([-2 28])
text(-1.75,1.7,['r = ' sprintf('%2.2f',rall)])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 2 -2 2])
title('All fishes')
% stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_',harv,'_SAUP_comp_types_temp_Stock_LELC_ppt.png'])
