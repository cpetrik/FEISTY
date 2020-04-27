% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, same k's
% Production & fisheries yield together on subplot

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
efile = 'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';

ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',efile,'/full_runs/'];

harv = 'All_fish03';

fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
epath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];

%% Yield
load([epath 'Hist_Fore_All_fish03_ensem6_mid_samek_ts_yield.mat'])
test=find(y>1950);
yid=test(1);

% types - diff
figure(1)
subplot(2,2,2)
plot(y,(rmF)*1e-6,'r','LineWidth',2); hold on;
plot(y,(rmP)*1e-6,'b','LineWidth',2); hold on;
plot(y,(rmD)*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
fill(X,(Rf)*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rp)*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rd)*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
text(y(yid),6.5,'B')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
legend('F','P','D')
legend('location','southwest')

% adults - diff
figure(2)
subplot(2,2,2)
plot(y,(rmF)*1e-6,'r','LineWidth',2); hold on;
plot(y,(rmLP)*1e-6,'b','LineWidth',2); hold on;
plot(y,(rmLD)*1e-6,'color',[0 0.6 0],'LineWidth',2); hold on;
fill(X,(Rf)*1e-6,'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rlp)*1e-6,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,(Rld)*1e-6,'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
xlim([y(yid) y(end)])
%ylim([-0.2 0.3])
text(y(yid),6.5,'B')
xlabel('Year')
ylabel('Yield (MT km^-^2 y^-^1) relative to 1951')
legend('F','P','D')
legend('location','southwest')

%% Production
load([epath2 'Hist_Fore_All_fish03_ensem6_mid_samek_ts_prod.mat'])
test=find(y>1950);
yid=test(1);

figure(1)
subplot(2,2,1)
fill(X,100*(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*(mpF),'r','LineWidth',2); hold on;
plot(y,100*(mpP),'b','LineWidth',2); hold on;
plot(y,100*(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
ylim([-35 15])
text(y(yid),0.18,'A')
xlabel('Year')
ylabel('% \Delta in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_yield_prod_types_all_ensem_mid6_samek_cone_1std_yr_2plot.png'])

figure(2)
subplot(2,2,1)
fill(X,100*(Vf),'r','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vp),'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*(Vd),'g','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*(mpF),'r','LineWidth',2); hold on;
plot(y,100*(mpP),'b','LineWidth',2); hold on;
plot(y,100*(mpD),'color',[0 0.6 0],'LineWidth',2); hold on;
xlim([y(yid) y(end)])
ylim([-35 15])
text(y(yid),0.18,'A')
xlabel('Year')
ylabel('% \Delta in production relative to 1951')
print('-dpng',[ppath 'Hist_Fore_',harv,'_biom_prod_types_adults_ensem_mid6_samek_cone_1std_yr_2plot.png'])





