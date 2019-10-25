% Time series of NPP, Z, and Det
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

%% means
mtemp = nanmean(ptemp);
mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_loss = nanmean(mz_loss);
mlz_loss = nanmean(lz_loss);
mL = nanmean(L);

%% change to integrated total each year
mz = mz_loss + lz_loss;
ZD = nanmean(mz ./ det);
lZD = nanmean(log10(mz ./ det)); 

%change to integrated total each year
% DO I NEED AREA INTEGRATED??? YES
area = AREA_OCN(ID);
area_mat = repmat(area,1,length(y));

tZD = nansum(mz.*area_mat) ./ nansum(det.*area_mat);
tlZD = log10( nansum(mz.*area_mat) ./ nansum(det.*area_mat) ); 

tPZ = nansum(mz.*area_mat) ./ nansum((mz.*area_mat)+(npp.*area_mat));
tPD = nansum(det.*area_mat) ./ nansum((det.*area_mat)+(npp.*area_mat));

tNPP = nansum(npp.*area_mat);
tZ = nansum(mz.*area_mat);
tD = nansum(det.*area_mat);

%% figure info
axesPosition = [110 40 200 200];  %# Axes position, in pixels
yWidth = 30;                      %# y axes spacing, in pixels
xLimit = [1895 2095];         %# Range of x values
xOffset = -yWidth*diff(xLimit)/axesPosition(3);

%% temp, npp, zoopL, det
%npp e-15 = Pg
%zoop & det e-12 = Tg

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cmap_ppt(5,:),...
          'XLim',xLimit,'YLim',[1.29 1.39],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(1,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[23 74],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[10 20],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tNPP(8:47)*1e-15,'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZ(8:47)*1e-12,'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),tD(8:47)*1e-12,'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_tot_npp_Zl_Det_temp.png'])

%% temp, npp, zoopL
%npp e-15 = Pg
%zoop & det e-12 = Tg

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cmap_ppt(5,:),...
          'XLim',xLimit,'YLim',[1.29 1.39],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(1,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[63 73],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[10 20],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tNPP(8:47)*1e-15,'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZ(8:47)*1e-12,'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_tot_npp_Zl_temp.png'])

%% temp, npp, det
%npp e-15 = Pg
%zoop & det e-12 = Tg

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cmap_ppt(5,:),...
          'XLim',xLimit,'YLim',[1.29 1.39],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(2,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[20 30],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14.5 18],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tNPP(8:47)*1e-15,'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tD(8:47)*1e-12,'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_tot_npp_Det_temp.png'])

%% npp, zoopL, det
%npp e-15 = Pg
%zoop & det e-12 = Tg

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cmap_ppt(5,:),...
          'XLim',xLimit,'YLim',[1.29 1.39],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(1,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[63 73],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(2,:),...
          'XLim',xLimit+[2*xOffset 0],'YLim',[20 30],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tNPP(8:47)*1e-15,'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZ(8:47)*1e-12,'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),tD(8:47)*1e-12,'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_tot_npp_Zl_Det.png'])

%% fractions zoopL/det, det/npp, zoopL/npp

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cm{6},...
          'XLim',xLimit,'YLim',[2.5 3],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(1,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[4.7 5.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(2,:),...
          'XLim',xLimit+[2*xOffset 0],'YLim',[1.7 2.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tPZ(8:47)*100,'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),tPD(8:47)*100,'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_frac_NZ_ND_ZlDet.png'])
