% Time series of P vs D and Z vs Det
% Area-integrated totals

clear all
close all

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

% GCM grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

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
load('cmap_ppt_angles.mat')
cmap3=cmap_ppt([1,3,5,4],:);

cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

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
% load(['Means_Historic_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
%     'mp_mean','md_mean','lp_mean','ld_mean','sf_mean','mf_mean');

HM = mf_mean+mp_mean+md_mean;
HL = lp_mean+ld_mean;
HF = sf_mean+mf_mean;
HP = sp_mean+mp_mean+lp_mean;
HD = sd_mean+md_mean+ld_mean;

clear sp_mean sd_mean mp_mean md_mean lp_mean ld_mean sf_mean mf_mean

%% RCP 8.5
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
   'mp_mean','md_mean','lp_mean','ld_mean','sf_mean','mf_mean');
% load(['Means_fore_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
%     'mp_mean','md_mean','lp_mean','ld_mean','sf_mean','mf_mean');

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
FD = nanmean(F ./ (F+D));
PelD = nanmean((P+F) ./ (P+F+D));
LM = nanmean(L ./ (L+M));

%change to integrated total each year
% DO I NEED AREA INTEGRATED??? YES
area = AREA_OCN(ID);
area_mat = repmat(area,1,length(y));

tZD = nansum(mz.*area_mat) ./ nansum(det.*area_mat);
tlZD = log10( nansum(mz.*area_mat) ./ nansum(det.*area_mat) ); 
tPD = nansum(P.*area_mat) ./ nansum((P.*area_mat)+(D.*area_mat));
tPF = nansum(P.*area_mat) ./ nansum((P.*area_mat)+(F.*area_mat));
tFD = nansum(F.*area_mat) ./ nansum((F.*area_mat)+(D.*area_mat));
tFP = nansum(F.*area_mat) ./ nansum((F.*area_mat)+(P.*area_mat));
tPelD = nansum((P.*area_mat)+(F.*area_mat)) ./ ...
    nansum((P.*area_mat)+(F.*area_mat)+(D.*area_mat));
tLM = nansum(L.*area_mat) ./ nansum((L.*area_mat)+(M.*area_mat));

%% means
mtemp = nanmean(ptemp);
mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_loss = nanmean(mz_loss);
mlz_loss = nanmean(lz_loss);
mL = nanmean(L);

%% figure info
axesPosition = [110 40 200 200];  %# Axes position, in pixels
yWidth = 30;                      %# y axes spacing, in pixels
xLimit = [1895 2095];         %# Range of x values
xOffset = -yWidth*diff(xLimit)/axesPosition(3);

%% Total global then divide; log10 ZlDet
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.3 0.8],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPelD(8:47),'color',[0 0.5 0.75],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_logZlDet_temp_all.png'])

%% no Pel:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.3 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_logZlDet_temp_noPelD.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
% h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
%           'Color','none','XColor','k','YColor','k',...
%           'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
%           'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
%line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_size_logZlDet.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_size_logZlDet_temp.png'])

%% Just Pel:Dem
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.77 0.79],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tPelD(8:47),'color',[0 0.5 0.75],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_PelDem_logZlDet_temp.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.77 0.79],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tPelD(8:47),'color',[0 0.5 0.75],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_PelDem_logZlDet.png'])

%% just F:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.65 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_logZlDet_temp.png'])

%% just F:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.65 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[0.4 0.45],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tlZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_logZlDet.png'])

%% Total global then divide; raw Zl:Det

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.8],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPelD(8:47),'color',[0 0.5 0.75],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_ZlDet_temp.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_ZlDet_temp_noPelDem.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
% h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
%           'Color','none','XColor','k','YColor','k',...
%           'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
%           'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
%line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_size_ZlDet.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tLM(8:47),'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPF(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPD(8:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_size_ZlDet_temp.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.77 0.79],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tPelD(8:47),'color',[0 0.5 0.75],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_PelDem_ZlDet_temp.png'])

%%
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.77 0.79],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tPelD(8:47),'color',[0 0.5 0.75],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_PelDem_ZlDet.png'])

%% F:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.65 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_ZlDet_temp.png'])

%% F:D vs Z:Det
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.65 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_ZlDet.png'])

%% F:P
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.54 0.64],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tFP(8:47),'color','r','Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FP_ZlDet_temp.png'])

%% Fake lines of expected trends
fakeLM = tLM(8)+ linspace(0.1,0,40);
fakePF = linspace(0.45,0.5,40);
fakePD = linspace(0.5,0.65,40);

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.35 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[2.5 2.9],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
% h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
%           'Color','none','XColor','k','YColor','k',...
%           'XLim',xLimit+[2*xOffset 0],'YLim',[14 18],...
%           'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),fakeLM,'color',cmap_ppt(2,:),'Linewidth',2,'LineStyle','--','Parent',h1); hold on;
line(y(8:47),fakePF,'color',cmap_ppt(1,:),'Linewidth',2,'LineStyle','--','Parent',h1); hold on;
line(y(8:47),fakePD,'color',cmap_ppt(5,:),'Linewidth',2,'LineStyle','--','Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_size_ZlDet_Fake.png'])


