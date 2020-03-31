% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Production & biomass together on subplot

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

%% time
load([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod.mat'],'X','y')
test=find(y>1950);
yid=test(1);

%% Biomass
load([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_biomass.mat'])

figure(1)
subplot(2,2,1)
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

subplot(2,2,3)
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
ylabel('% \Delta in biomass relative to 1951')

%% Production
load([epath 'Hist_Fore_All_fish03_ensem6_mid_temp3_ts_prod.mat'])

subplot(2,2,2)
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

subplot(2,2,4)
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
ylabel('% \Delta in production relative to 1951')

print('-dpng',[ppath 'Hist_Fore_',harv,'_biom_prod_types_ensem_mid6_temp3_cone_1std_yr_4plot.png'])





