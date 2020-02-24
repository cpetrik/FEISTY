% Time series of P vs D and Z vs Det
% Area-integrated totals
% Z prod instead off loss

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

% GCM grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

% Zoop and det and npp 
load([gpath 'cobalt_det_temp_zoop_npp_means.mat']);

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

%% Original parameters
lpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
%fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/',cfile '/'];
fpath = [lpath,cfile '/'];
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat']);

%% Ensemble parameter sets
% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = [lpath,...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);

%% ts
%In original saved file
% y1 = 1860+(1/12):(1/12):2005;
% y2 = 2005+(1/12):(1/12):2100;
% y = [y1 y2];

% SOMETHING WRONG WITH MOVING MEAN OF #28 AROUND 1949
% FIX TS OF SMALL FISH AT T=1080
%tF[1081,28]=1.94e+35
%tP[1081,28]=1.94e+35
%tD[1081,28]=1.94e+35
%tA[1081,28]=5.83e+35
hTsF(28,1081) = (hTsF(28,1080) + hTsF(28,1082))/2;
hTsP(28,1081) = (hTsP(28,1080) + hTsP(28,1082))/2;
hTsD(28,1081) = (hTsD(28,1080) + hTsD(28,1082))/2;

HF = hTsF + hTmF;
HP = hTsP + hTmP + hTlP;
HD = hTsD + hTmD + hTlD;
HA = HF + HP + HD;

FF = fTsF + fTmF;
FP = fTsP + fTmP + fTlP;
FD = fTsD + fTmD + fTlD;
FA = FF + FP + FD;

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];

npp = [hmnpp fmnpp];
det = [hmdet fmdet];
mz_prod = [hmmz_prod fmmz_prod];
lz_prod = [hmlz_prod fmlz_prod];
ptemp = [ptemp_5yr_hist ptemp_5yr_fore] - 273;

%% ratios and fractions
mz = mz_prod + lz_prod;

%change to integrated total each year
% DO I NEED AREA INTEGRATED??? YES
area = AREA_OCN(ID);
area_mat = repmat(area,1,length(y));

tZD = nansum(mz.*area_mat) ./ nansum(det.*area_mat);
tPD = nansum(P.*area_mat) ./ nansum((P.*area_mat)+(D.*area_mat));
tFD = nansum(F.*area_mat) ./ nansum((F.*area_mat)+(D.*area_mat));
tPelD = nansum((P.*area_mat)+(F.*area_mat)) ./ ...
    nansum((P.*area_mat)+(F.*area_mat)+(D.*area_mat));

%% CONE OF UNCERTAINTY
% means
mtemp = nanmean(ptemp);
mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_prod = nanmean(mz_prod);
mlz_prod = nanmean(lz_prod);
mL = nanmean(L);
mHTL = mean([eTE_HTL; q(5,:)]);
mATL = mean([eTE_ATL; q(6,:)]);
sHTL = std([eTE_HTL; q(5,:)]);
sATL = std([eTE_ATL; q(6,:)]);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%+/- 1 stdev
Sa=[mATL+sATL fliplr(mATL-sATL)]; 
Sh=[mHTL+sATL fliplr(mHTL-sATL)]; 

%% figure info
axesPosition = [110 40 200 200];  %# Axes position, in pixels
yWidth = 30;                      %# y axes spacing, in pixels
xLimit = [1895 2095];         %# Range of x values
xOffset = -yWidth*diff(xLimit)/axesPosition(3);

%% PvD with ZprodDet, no temp 
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.55 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[4.7 5.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tPD(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
print('-dpng',[ppath 'Hist_Fore_',harv,'_PDfrac_ZpDet.png'])

%% PvD with ZprodDet, with temp
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.55 0.65],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[4.7 5.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14.5 17.5],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tPD(8:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
print('-dpng',[ppath 'Hist_Fore_',harv,'_PDfrac_ZpDet_temp.png'])

%% Just Pel:Dem
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.77 0.79],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[4.7 5.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14.5 17.5],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tPelD(8:47),'color',[0.5 0.5 0.5],'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
print('-dpng',[ppath 'Hist_Fore_',harv,'_PelDfrac_ZpDet_temp.png'])

%% just F:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.65 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[4.7 5.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[14.5 17.5],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(8:47),mtemp(8:47),'color','k','Linewidth',2,'Parent',h3); hold on;
print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_logZlDet_temp.png'])

%% just F:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[0.65 0.7],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[4.7 5.2],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(8:47),tFD(8:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(8:47),tZD(8:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
%print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_logZlDet.png'])



