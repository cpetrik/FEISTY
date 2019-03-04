% Time series of PvsD and ZvsDet

clear all
close all

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% Zoop and det and npp 
load([gpath 'cobalt_det_temp_zoop_npp_means.mat']);

% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

cmYOR=cbrewer('seq','YlOrRd',50);
cmRP=cbrewer('seq','RdPu',50);
cmPR=cbrewer('seq','PuRd',50);


%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mzloss_5yr_hist = mzloss_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_5yr_hist = lzloss_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_hist = det_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_hist = npp_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

hmmz_loss = mzloss_5yr_hist;
hmlz_loss = lzloss_5yr_hist;
hmdet = det_5yr_hist;
hmnpp = npp_5yr_hist;

mzloss_5yr_fore = mzloss_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_5yr_fore = lzloss_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_fore = det_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_fore = npp_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

fmmz_loss = mzloss_5yr_fore;
fmlz_loss = lzloss_5yr_fore;
fmdet = det_5yr_fore;
fmnpp = npp_5yr_fore;

%% Hindcast
load([fpath 'Means_Historic_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
    'mp_mean','md_mean','lp_mean','ld_mean','sf_mean','mf_mean');

HM = mf_mean+mp_mean+md_mean;
HL = lp_mean+ld_mean;
HF = sf_mean+mf_mean;
HP = sp_mean+mp_mean+lp_mean;
HD = sd_mean+md_mean+ld_mean;

clear sp_mean sd_mean mp_mean md_mean lp_mean ld_mean sf_mean mf_mean

%% RCP 8.5
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
    'mp_mean','md_mean','lp_mean','ld_mean','sf_mean','mf_mean');

FM = mf_mean+mp_mean+md_mean;
FL = lp_mean+ld_mean;
FF = sf_mean+mf_mean;
FP = sp_mean+mp_mean+lp_mean;
FD = sd_mean+md_mean+ld_mean;

clear sp_mean sd_mean mp_mean md_mean lp_mean ld_mean sf_mean mf_mean

%% ts
y1 = 1860+(1/12):5:2005;
y2 = 2005+(1/12):5:2100;
y = [y1 y2];

M = [HM FM];
L = [HL FL];
F = [HF FF];
P = [HP FP];
D = [HD FD];
npp = [hmnpp fmnpp];
det = [hmdet fmdet];
mz_loss = [hmmz_loss fmmz_loss];
lz_loss = [hmlz_loss fmlz_loss];
ptemp = [ptemp_5yr_hist ptemp_5yr_fore] - 273;

%% ratios and fractions
mz = mz_loss + lz_loss;
ZD = nanmean(mz ./ det);
lZD = nanmean(log10(mz ./ det));

PD = nanmean(P ./ (P+D));
PF = nanmean(P ./ (P+F));
PelD = nanmean((P+F) ./ (P+F+D));
LM = nanmean(L ./ (L+M));

%% means
mtemp = nanmean(ptemp);
mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_loss = nanmean(mz_loss);
mlz_loss = nanmean(lz_loss);
mL = nanmean(L);


%%
figure(1)
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(PelD),'b','Linewidth',2); hold on;
legend('P/(P+F)','P/(P+D)','Pelag/(Pelag+Dem)')
legend('location','eastoutside')
xlim([1895 2095])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_type_v2.png'])

figure(2)
plot(y,(LM),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(PelD),'b','Linewidth',2); hold on;
legend('L/(L+M)','P/(P+F)','P/(P+D)','Pelag/(Pelag+Dem)')
legend('location','eastoutside')
xlim([1895 2095])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size_v2.png'])

%% no legends
figure(3)
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(PelD),'b','Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_type_v2_noleg.png'])

figure(4)
plot(y,(LM),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(PF),'r','Linewidth',2); hold on;
plot(y,(PD),'k','Linewidth',2); hold on;
plot(y,(PelD),'b','Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size_v2_noleg.png'])

%%
figure(5)
plot(mtemp,(LM),'.','color',[0.5 0.5 0.5],'MarkerSize',20); hold on;
plot(mtemp,(PF),'.r','MarkerSize',20); hold on;
plot(mtemp,(PD),'.k','MarkerSize',20); hold on;
plot(mtemp,(PelD),'.b','MarkerSize',20); hold on;
legend('L/(L+M)','P/(P+F)','P/(P+D)','Pelag/(Pelag+Dem)')
legend('location','east')
xlabel('Mean Temperature')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_temp.png'])

%%
figure(6)
plot(ZD,(LM),'.','color',[0.5 0.5 0.5],'MarkerSize',20); hold on;
plot(ZD,(PF),'.r','MarkerSize',20); hold on;
plot(ZD,(PD),'.k','MarkerSize',20); hold on;
plot(ZD,(PelD),'.b','MarkerSize',20); hold on;
legend('L/(L+M)','P/(P+F)','P/(P+D)','Pelag/(Pelag+Dem)')
legend('location','east')
xlabel('Zoop : Det')
% ylabel('Fractions')
% title(['Fished'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_ZlDet.png'])

%%
figure(7)
subplot(2,2,1)
plot(mtemp,(LM),'.','color',[0.5 0.5 0.5],'MarkerSize',20); hold on;
plot(mtemp,(PF),'.r','MarkerSize',20); hold on;
plot(mtemp,(PD),'.k','MarkerSize',20); hold on;
plot(mtemp,(PelD),'.b','MarkerSize',20); hold on;
ylim([0.1 0.9])
xlabel('Mean Temperature')

subplot(2,2,2)
plot(ZD,(LM),'.','color',[0.5 0.5 0.5],'MarkerSize',20); hold on;
plot(ZD,(PF),'.r','MarkerSize',20); hold on;
plot(ZD,(PD),'.k','MarkerSize',20); hold on;
plot(ZD,(PelD),'.b','MarkerSize',20); hold on;
ylim([0.1 0.9])
xlabel('Zoop : Det')
print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_temp_ZlDet.png'])

%% Mult y axes
f8=figure(8);
lc = [0 0.5 0.5];
rc = [0.5 0.5 0.5];
set(f8,'defaultAxesColorOrder',[lc; rc]);
yyaxis left
plot(y,(LM),'Linewidth',2); hold on;
plot(y,(PF),'Linewidth',2); hold on;
plot(y,(PD),'Linewidth',2); hold on;
plot(y,(PelD),'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fractions')

yyaxis right
plot(y,mtemp,'Linewidth',2); hold on;
plot(y,ZD,'Linewidth',2); hold on;
xlim([1895 2095])
ylabel('Temperature and Z:D')
%print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size_v2_noleg.png'])

%%
f9=figure(9);
set(f9,'defaultAxesColorOrder',[lc; rc]);
yyaxis left
plot(y,(LM),'Linewidth',2); hold on;
plot(y,(PF),'Linewidth',2); hold on;
plot(y,(PD),'Linewidth',2); hold on;
plot(y,(PelD),'Linewidth',2); hold on;
plot(y,lZD,'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fractions and log_1_0Z:D')

yyaxis right
plot(y,mtemp,'Linewidth',2); hold on;
xlim([1895 2095])
ylabel('Temperature')
%print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size_v2_noleg.png'])

%%
load('cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm={[1 0.5 0],...   %1 orange
    [0.5 0.5 0],... %2 tan/army
    [0 0.7 0],...   %3 g
    [0 1 1],...     %4 c
    [0 0 0.75],...  %5 b
    [0.5 0 1],...   %6 purple
    [1 0 1],...     %7 m
    [1 0 0],...     %8 r
    [0.5 0 0],...   %9 maroon
    [0.75 0.75 0.75],... %10 lt grey
    [0.5 0.5 0.5],...    %11 med grey
    [49/255 79/255 79/255],... %12 dk grey
    [0 0 0],...      %13 black
    [1 1 0],...      %14 yellow
    [127/255 255/255 0],... %15 lime green
    [0 0.5 0],...    %16 dk green
    [0/255 206/255 209/255],... %17 turq
    [0 0.5 0.75],...   %18 med blue
    [188/255 143/255 143/255],... %19 rosy brown
    [255/255 192/255 203/255],... %20 pink
    [255/255 160/255 122/255]}; %21 peach

f10=figure(10);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
plot(y,(LM),'color',cmap_ppt(1,:),'Linewidth',2); hold on;
plot(y,(PF),'color',cmap_ppt(3,:),'Linewidth',2); hold on;
plot(y,(PD),'color',cmap_ppt(5,:),'Linewidth',2); hold on;
plot(y,(PelD),'color',cmap_ppt(2,:),'Linewidth',2); hold on;
plot(y,lZD,'color',cm{6},'Linewidth',2); hold on;
plot(y,lZD,'color',[0.5 0 0.5],'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fractions and log_1_0Z:D')

ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'Xlim',[1895 2095],...
    'Ylim',[14 18]);
line(y,mtemp,'color','k','Linewidth',2,'Parent',ax2); hold on;
ylabel('Temperature')
%print('-dpng',[ppath 'Hist_Fore_',harv,'_fracs_size_v2_noleg.png'])

%%
figure(11);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
plot(y,(LM),'color',cmap_ppt(1,:),'Linewidth',2); hold on;
plot(y,(PF),'color',cmap_ppt(3,:),'Linewidth',2); hold on;
plot(y,(PD),'color',cmap_ppt(5,:),'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fractions')

ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'Xlim',[1895 2095]);
line(y,mtemp,'color','k','Linewidth',2,'Parent',ax2); hold on;
ylabel('Temperature')
%%
figure(12);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
plot(y,(PelD),'color',[0 0.5 0.75],'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fraction all pelagics (F+P)')

ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'Xlim',[1895 2095]);
line(y,mtemp,'color','k','Linewidth',2,'Parent',ax2); hold on;
ylabel('Temperature')
%%
figure(13);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
plot(y,(PelD),'color',[0 0.5 0.75],'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fraction')

ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'Xlim',[1895 2095]);
line(y,lZD,'color',cm{6},'Linewidth',2,'Parent',ax2); hold on;
ylabel('log_1_0Z:D')
%%
figure(14);
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
plot(y,(LM),'color',cmap_ppt(1,:),'Linewidth',2); hold on;
plot(y,(PF),'color',cmap_ppt(3,:),'Linewidth',2); hold on;
plot(y,(PD),'color',cmap_ppt(5,:),'Linewidth',2); hold on;
xlim([1895 2095])
xlabel('Year')
ylabel('Fractions and log_1_0Z:D')

ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'Xlim',[1895 2095]);
line(y,lZD,'color',cm{6},'Linewidth',2,'Parent',ax2); hold on;
ylabel('log_1_0Z:D')

%%
figure(15);
plot(y,(LM),'color',cmap_ppt(1,:),'Linewidth',2); hold on;
plot(y,(PF),'color',cmap_ppt(3,:),'Linewidth',2); hold on;
plot(y,(PD),'color',cmap_ppt(5,:),'Linewidth',2); hold on;
plotyyy(y,(PelD),y,lZD,y,mtemp,{'Fractions','log_1_0Z:D','Temperature'});
xlim([1895 2095])
xlabel('Year')
