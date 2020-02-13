% Time series by type and biome
% Use historic biomes into future
% Visualize time series output of FEISTY
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files

clear all
close all

bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);
load([bpath 'COBALT_biomes_50yr_means_5100.mat'],'biome_hist','biome_fore');

gb_hist=biome_hist(grid(:,1));
gb_fore=biome_fore(grid(:,1));
hLC = (gb_hist==1);
hCS = (gb_hist==2);
hSS = (gb_hist==3);
fLC = (gb_fore==1);
fCS = (gb_fore==2);
fSS = (gb_fore==3);

% Colors
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmapT = cmR(35,:);
cmapT(2,:) = cmB(35,:);
cmapT(3,:) = cmG(35,:);
cmapT(4,:)=[0 0 0];

load('cmap_ppt_angles.mat')
cmap3(1,:)=cmap_ppt(1,:);
cmap3(2,:)=cmap_ppt(3,:);
cmap3(3,:)=cmap_ppt(5,:);

%% FEISTY output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_fore_',harv,'_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'b_mean','lp_mean','ld_mean');
Ft5.F = sf_mean+mf_mean;
Ft5.P = sp_mean+mp_mean+lp_mean;
Ft5.D = sd_mean+md_mean+ld_mean;
Ft5.B = b_mean;
Ft5.A = Ft5.F + Ft5.P + Ft5.D;
Ft5.M = mf_mean+mp_mean+md_mean;
Ft5.L = lp_mean+ld_mean;
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean b_mean


load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'b_mean','lp_mean','ld_mean');
Ht5.F = sf_mean+mf_mean;
Ht5.P = sp_mean+mp_mean+lp_mean;
Ht5.D = sd_mean+md_mean+ld_mean;
Ht5.B = b_mean;
Ht5.A = Ht5.F + Ht5.P + Ht5.D;
Ht5.M = mf_mean+mp_mean+md_mean;
Ht5.L = lp_mean+ld_mean;
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean b_mean

% save hist and fore together
save([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat'],'-append')

%% ts
y1 = 1860+(1/12):5:2005;
y2 = 2005+(1/12):5:2100;
y = [y1 y2];

mFF = [nanmean(Ht5.F) nanmean(Ft5.F)];
mFP = [nanmean(Ht5.P) nanmean(Ft5.P)];
mFD = [nanmean(Ht5.D) nanmean(Ft5.D)];
mFA = [nanmean(Ht5.A) nanmean(Ft5.A)];

mLCF = [nanmean(Ht5.F(hLC,:)) nanmean(Ft5.F(fLC,:))];
mLCP = [nanmean(Ht5.P(hLC,:)) nanmean(Ft5.P(fLC,:))];
mLCD = [nanmean(Ht5.D(hLC,:)) nanmean(Ft5.D(fLC,:))];
mLCA = [nanmean(Ht5.A(hLC,:)) nanmean(Ft5.A(fLC,:))];
% mLCP needs smoothing, discontinuity at 2000
mmLCP = movmean(mLCP,10);
% splice together
mLCP2 = mLCP;
mLCP2(25:35) = mmLCP(25:35);

mCSF = [nanmean(Ht5.F(hCS,:)) nanmean(Ft5.F(fCS,:))];
mCSP = [nanmean(Ht5.P(hCS,:)) nanmean(Ft5.P(fCS,:))];
mCSD = [nanmean(Ht5.D(hCS,:)) nanmean(Ft5.D(fCS,:))];
mCSA = [nanmean(Ht5.A(hCS,:)) nanmean(Ft5.A(fCS,:))];

mSSF = [nanmean(Ht5.F(hSS,:)) nanmean(Ft5.F(fSS,:))];
mSSP = [nanmean(Ht5.P(hSS,:)) nanmean(Ft5.P(fSS,:))];
mSSD = [nanmean(Ht5.D(hSS,:)) nanmean(Ft5.D(fSS,:))];
mSSA = [nanmean(Ht5.A(hSS,:)) nanmean(Ft5.A(fSS,:))];

%%
figure(1)
plot(y,log10(mFF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mFP),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mFD),'color',cmapT(3,:),'Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-0.15 0.3])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ts5.png'])

%%
figure(2)
subplot(3,3,1)
plot(y,log10(mLCF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mLCP2),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mLCD),'color',cmapT(3,:),'Linewidth',2); hold on;
% legend('F','P','D')
% legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-1 0.5])
title('LC')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,2)
plot(y,log10(mCSF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSP),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mCSD),'color',cmapT(3,:),'Linewidth',2); hold on;
% legend('F','P','D')
% legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-1 0.5])
title('ECCS')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,3)
plot(y,log10(mSSF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mSSP),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSD),'color',cmapT(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-1 0.5])
title('ECSS')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,7)
plot(y,log10(mLCF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mLCP2),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mLCD),'color',cmapT(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-2 0.25])
title('LC')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,8)
plot(y,log10(mCSF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSP),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mCSD),'color',cmapT(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
%ylim([-1 0.5])
title('ECCS')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,9)
plot(y,log10(mSSF),'color',cmapT(1,:),'Linewidth',2); hold on;
plot(y,log10(mSSP),'color',cmapT(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSD),'color',cmapT(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
%ylim([-1 0.5])
title('ECSS')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ts5_biomes.png'])

%%
figure(3)
subplot(3,3,1)
plot(y,log10(mLCF),'color',cmap3(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSF),'color',cmap3(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSF),'color',cmap3(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-1 0.5])
title('F')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,2)
plot(y,log10(mLCP2),'color',cmap3(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSP),'color',cmap3(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSP),'color',cmap3(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-1 0.5])
title('P')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,3)
plot(y,log10(mLCD),'color',cmap3(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSD),'color',cmap3(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSD),'color',cmap3(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-1 0.5])
title('D')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,7)
plot(y,log10(mLCF),'color',cmap3(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSF),'color',cmap3(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSF),'color',cmap3(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-0.25 0.5])
title('F')
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,8)
plot(y,log10(mLCP2),'color',cmap3(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSP),'color',cmap3(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSP),'color',cmap3(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
ylim([-2 0.5])
title('P')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

subplot(3,3,9)
plot(y,log10(mLCD),'color',cmap3(1,:),'Linewidth',2); hold on;
plot(y,log10(mCSD),'color',cmap3(2,:),'Linewidth',2); hold on;
plot(y,log10(mSSD),'color',cmap3(3,:),'Linewidth',2); hold on;
xlim([y(1) y(end)])
%ylim([-1 0.5])
title('D')
xlabel('Year')
set(gca,'XTick',[1950 2050],'XTickLabel',[1950 2050])

print('-dpng',[ppath 'Hist_Fore_',harv,'_all_biomes_ts5_types.png'])





