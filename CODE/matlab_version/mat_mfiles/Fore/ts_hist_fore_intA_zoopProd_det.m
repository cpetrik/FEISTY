% Time series of NPP, Z, and Det
% Area-integrated totals

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
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
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

cmYOR=cbrewer('seq','YlOrRd',50,'PCHIP');
cmRP=cbrewer('seq','RdPu',50,'PCHIP');
cmPR=cbrewer('seq','PuRd',50,'PCHIP');
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
mzprod_5yr_hist = mzprod_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_5yr_hist = lzprod_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_hist = det_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_hist = npp_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

hmmz_prod = mzprod_5yr_hist;
hmlz_prod = lzprod_5yr_hist;
hmdet = det_5yr_hist;
hmnpp = npp_5yr_hist;

mzprod_5yr_fore = mzprod_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_5yr_fore = lzprod_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_fore = det_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_fore = npp_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

fmmz_prod = mzprod_5yr_fore;
fmlz_prod = lzprod_5yr_fore;
fmdet = det_5yr_fore;
fmnpp = npp_5yr_fore;

%% ts
y1 = 1860+(1/12):5:2005;
y2 = 2005+(1/12):5:2100;
y = [y1 y2];

npp = [hmnpp fmnpp];
det = [hmdet fmdet];
mz_prod = [hmmz_prod fmmz_prod];
lz_prod = [hmlz_prod fmlz_prod];
ptemp = [ptemp_5yr_hist ptemp_5yr_fore] - 273;

%% means
mtemp = nanmean(ptemp);
mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_prod = nanmean(mz_prod);
mlz_prod = nanmean(lz_prod);

%% change to integrated total each year
mz = mz_prod + lz_prod;
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
xLimit = [1950 2095];         %# Range of x values
xOffset = -yWidth*diff(xLimit)/axesPosition(3);

%% fractions zoopL/det, det/npp, zoopL/npp

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cm{6},...
          'XLim',xLimit,'YLim',[4.6 5.3],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(1,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[0.084e2 0.091e2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(2,:),...
          'XLim',xLimit+[2*xOffset 0],'YLim',[0.015e2 0.022e2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(20:47),tZD(20:47),'color',cm{6},'Linewidth',2,'Parent',h1); hold on;
line(y(20:47),tPZ(20:47)*100,'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h2); hold on;
line(y(20:47),tPD(20:47)*100,'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
text(1937,2.25,'Z:D','color',cm{6});
text(1913,2.25,'Z:P','color',cmap_ppt(1,:));
text(1889,2.25,'D:P','color',cmap_ppt(2,:));
print('-dpng',[ppath 'Hist_Fore_frac_ZpN_DN_ZpDet.png'])

%% zoop, det, fraction zoopP/det

figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cm{6},...
          'XLim',xLimit,'YLim',[4.7 5.2],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(1,:),...
          'XLim',xLimit+[xOffset 0],'YLim',[1.0 1.5],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cmap_ppt(2,:),...
          'XLim',xLimit+[2*xOffset 0],'YLim',[2.3 2.8],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');
%ylabel(h3,'values');

line(y(20:47),tZD(20:47),'color',cm{6},'Linewidth',2,'Parent',h1); hold on;
line(y(20:47),tZ(20:47)*1e-14,'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h2); hold on;
line(y(20:47),tD(20:47)*1e-13,'color',cmap_ppt(2,:),'Linewidth',2,'Parent',h3); hold on;
% xlim([1895 2095])
%xlabel('Year')
text(1937,2.83,'Z:D','color',cm{6});
text(1917,2.83,'Z','color',cmap_ppt(1,:));
text(1897,2.83,'D','color',cmap_ppt(2,:));
print('-dpng',[ppath 'Hist_Fore_Zp_D_ZpDet.png'])
