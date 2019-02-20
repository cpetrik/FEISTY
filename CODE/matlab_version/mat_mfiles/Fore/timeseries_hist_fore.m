% Visualize time series output of FEISTY
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'ForeFish');
load([fpath 'Means_fore_pristine_' cfile '.mat'],'ForePris');
load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],'HistFish');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% ts
y1 = 1860+(1/12):(1/12):2005;
y2 = 2005+(1/12):(1/12):2100;
y = [y1 y2];

HFF = HistFish(1,:)+HistFish(4,:);
HFP = HistFish(2,:)+HistFish(5,:)+HistFish(7,:);
HFD = HistFish(3,:)+HistFish(6,:)+HistFish(8,:);
HFA = HFF + HFP + HFD;

FPF = ForePris(1,:)+ForePris(4,:);
FPP = ForePris(2,:)+ForePris(5,:)+ForePris(7,:);
FPD = ForePris(3,:)+ForePris(6,:)+ForePris(8,:);
FPA = FPF + FPP + FPD;

FFF = ForeFish(1,:)+ForeFish(4,:);
FFP = ForeFish(2,:)+ForeFish(5,:)+ForeFish(7,:);
FFD = ForeFish(3,:)+ForeFish(6,:)+ForeFish(8,:);
FFA = FFF + FFP + FFD;

%%
figure(1)
plot(y1,log10(HFF),'r','Linewidth',2); hold on;
plot(y1,log10(HFP),'b','Linewidth',2); hold on;
plot(y1,log10(HFD),'k','Linewidth',2); hold on;
plot(y2,log10(FPF),':r','Linewidth',2); hold on;
plot(y2,log10(FPP),':b','Linewidth',2); hold on;
plot(y2,log10(FPD),':k','Linewidth',2); hold on;
plot(y2,log10(FFF),'r','Linewidth',2); hold on;
plot(y2,log10(FFP),'b','Linewidth',2); hold on;
plot(y2,log10(FFD),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.15 0.3])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
%title(['Historic fished'])
print('-dpng',[ppath 'Hist_Fore_pristine_',harv,'_all_types.png'])

%%
figure(2)
plot(y1,log10(HFA),'k','LineWidth',2); hold on;
plot(y2,log10(FPA),':k','LineWidth',2); hold on;
plot(y2,log10(FFA),'k','LineWidth',2); hold on;
ylim([0.46 0.62])
xlim([y(1) y(end)])
xlabel('Year')
ylabel('All fish mean biomass (g/m^2)')
%title('Forecast fished')
print('-dpng',[ppath 'Hist_Fore_pristine_',harv,'_ts_mbio.png'])

%%
FF = [HFF FFF];
FP = [HFP FFP];
FD = [HFD FFD];
FA = [HFA FFA];
figure(3)
plot(y,log10(FF),'r','Linewidth',2); hold on;
plot(y,log10(FP),'b','Linewidth',2); hold on;
plot(y,log10(FD),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.15 0.3])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types.png'])

%%
figure(4)
plot(y,log10(FA),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
ylim([0.46 0.62])
xlim([y(1) y(end)])
xlabel('Year')
ylabel('All fish mean biomass (g/m^2)')
title('Fished')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ts_mbio.png'])

test = movmean(FA,61);
hold on
plot(y,log10(test),'c','LineWidth',2);



