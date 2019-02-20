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
load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],'HistFish');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% ts
y1 = 1860+(1/12):(1/12):2005;
y2 = 2005+(1/12):(1/12):2100;
y = [y1 y2];

HF = HistFish(1,:)+HistFish(4,:);
HP = HistFish(2,:)+HistFish(5,:)+HistFish(7,:);
HD = HistFish(3,:)+HistFish(6,:)+HistFish(8,:);
HA = HF + HP + HD;
HM = sum(HistFish(4:6,:));
HL = sum(HistFish(7:8,:));
Hpel = HF + HP;

FF = ForeFish(1,:)+ForeFish(4,:);
FP = ForeFish(2,:)+ForeFish(5,:)+ForeFish(7,:);
FD = ForeFish(3,:)+ForeFish(6,:)+ForeFish(8,:);
FA = FF + FP + FD;
FM = sum(ForeFish(4:6,:));
FL = sum(ForeFish(7:8,:));
Fpel = FF + FP;

F = [HF FF];
P = [HP FP];
D = [HD FD];
A = [HA FA];
M = [HM FM];
L = [HL FL];
Pel = [Hpel Fpel];

PD = P ./ (P+D);
PF = P ./ (P+F);
LM = L ./ (L+M);
pd = Pel ./ (Pel+D);

%%
figure(1)
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(pd),'b','Linewidth',2); hold on;
legend('P/(P+F)','P/(P+D)','Pelag/(Pelag+Dem)')
legend('location','eastoutside')
xlim([1900 2100])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_type.png'])

figure(2)
plot(y,(LM),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(pd),'b','Linewidth',2); hold on;
legend('L/(L+M)','P/(P+F)','P/(P+D)','Pelag/(Pelag+Dem)')
legend('location','eastoutside')
xlim([1900 2100])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size.png'])

%% no legends
figure(3)
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(pd),'b','Linewidth',2); hold on;
xlim([1900 2100])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_type_noleg.png'])

figure(4)
plot(y,(LM),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(pd),'b','Linewidth',2); hold on;
xlim([1900 2100])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size_noleg.png'])





