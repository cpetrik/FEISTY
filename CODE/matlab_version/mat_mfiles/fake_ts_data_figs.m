% Fake Time series of ensembles
% Temp, O2, plankton

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
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
epath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

ppath = '/Users/cpetrik/Dropbox/ImptDocs/Funding/NOAA/ClimatePrograms/NCAR/schematic/';

%% colors
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

%% Hindcast
load([fpath 'Historic_ESM2M/Means_Historic_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
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
load([fpath 'Forecast_RCP85_ESM2M/Means_fore_',harv,'_' cfile '.mat'],'sp_mean','sd_mean',...
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
mz_prod = [hmmz_prod fmmz_prod];
lz_prod = [hmlz_prod fmlz_prod];
ptemp = [ptemp_5yr_hist ptemp_5yr_fore] - 273;

%% ratios and fractions
mz = mz_prod + lz_prod;
All = F + P + D;
ZD = nanmean(mz ./ det);
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
tPD = nansum(P.*area_mat) ./ nansum((P.*area_mat)+(D.*area_mat));
tPF = nansum(P.*area_mat) ./ nansum((P.*area_mat)+(F.*area_mat));
tFD = nansum(F.*area_mat) ./ nansum((F.*area_mat)+(D.*area_mat));
tFP = nansum(F.*area_mat) ./ nansum((F.*area_mat)+(P.*area_mat));
tPelD = nansum((P.*area_mat)+(F.*area_mat)) ./ ...
    nansum((P.*area_mat)+(F.*area_mat)+(D.*area_mat));
tLM = nansum(L.*area_mat) ./ nansum((L.*area_mat)+(M.*area_mat));

tP = nansum(P.*area_mat) ./ nansum(All.*area_mat);
tF = nansum(F.*area_mat) ./ nansum(All.*area_mat);
tPel = nansum((P.*area_mat)+(F.*area_mat)) ./ nansum(All.*area_mat);

%% means
mtemp = nanmean(ptemp);
mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_prod = nanmean(mz_prod);
mlz_prod = nanmean(lz_prod);
mL = nanmean(L);

%% diff from 1951
id = find(y>1950);
yid = id(1);
tFD = tFD - tFD(yid);
tPD = tPD - tPD(yid);
tPF = tPF - tPF(yid);
tZD = tZD - tZD(yid);
mtemp = mtemp - mtemp(yid);

%% fakes
time = y(19:47);
plank = rescale(tPD(19:47),0,1);
oxy = rescale(tPF(19:47),0,1);
temp1 = rescale(mtemp(19:47),0,1);
temp2 = rescale(tZD(19:47),0,1);
temp3 = rescale(tFD(19:47),0,1);

r1 = randn(length(time),10);
r2 = randn(length(time),10);
r3 = randn(length(time),10);

eplank = 0.25*r1' + plank;
eoxy = 0.25*r2' + oxy;
etemp1 = 0.25*r3' + temp1;
etemp2 = 0.25*r3' + temp2;
etemp3 = 0.25*r3' + temp3;

%% save

save([ppath 'ensemble_timeseries_data.mat'],'plank','eplank','oxy',...
    'eoxy','temp1','etemp1','time');

Tvec(1,:) = time;
Tvec(2,:) = temp1;
Tvec(3:12,:) = etemp1;
Tvec = Tvec';
Ttab = array2table(Tvec,'VariableNames',{'time','temp1','temp2','temp3','temp4',...
    'temp5','temp6','temp7','temp8','temp9','temp10','temp11'});
writetable(Ttab,[ppath 'temp_ts_data.csv'],'Delimiter',',')

Pvec(1,:) = time;
Pvec(2,:) = plank;
Pvec(3:12,:) = eplank;
Pvec = Pvec';
Ptab = array2table(Pvec,'VariableNames',{'time','plank1','plank2','plank3','plank4',...
    'plank5','plank6','plank7','plank8','plank9','plank10','plank11'});
writetable(Ptab,[ppath 'plank_ts_data.csv'],'Delimiter',',')

Ovec(1,:) = time;
Ovec(2,:) = oxy;
Ovec(3:12,:) = eoxy;
Ovec = Ovec';
Otab = array2table(Ovec,'VariableNames',{'time','oxy1','oxy2','oxy3','oxy4',...
    'oxy5','oxy6','oxy7','oxy8','oxy9','oxy10','oxy11'});
writetable(Otab,[ppath 'oxy_ts_data.csv'],'Delimiter',',')


%% figure info
axesPosition = [110 40 200 200];  %# Axes position, in pixels
yWidth = 30;                      %# y axes spacing, in pixels
xLimit = [1950 2095];         %# Range of x values
xOffset = -yWidth*diff(xLimit)/axesPosition(3);


%% PvD & PvF with ZprodDet, with temp
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[-0.06 0.02],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[-0.1 0.4],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[-0.2 1.8],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(19:47),tPD(19:47),'color',cmap_ppt(5,:),'Linewidth',2,'Parent',h1); hold on;
line(y(19:47),tPF(19:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
line(y(19:47),tZD(19:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(19:47),mtemp(19:47),'color','k','Linewidth',2,'Parent',h3); hold on;
%print('-dpng',[ppath 'Hist_Fore_',harv,'_PDfrac_PFfrac_ZpDet_temp_diff1951.png'])

%% just F:D
figure('Units','pixels','Position',[200 200 330 260]);
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
          'XLim',xLimit,'YLim',[-0.005 0.03],'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',[0.5 0 1.0],...
          'XLim',xLimit+[xOffset 0],'YLim',[-0.1 0.4],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'YLim',[-0.2 1.8],...
          'XTick',[],'XTickLabel',[],'NextPlot','add');
xlabel(h1,'Year');

line(y(19:47),tFD(19:47),'color',cmap_ppt(3,:),'Linewidth',2,'Parent',h1); hold on;
line(y(19:47),tZD(19:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
line(y(19:47),mtemp(19:47),'color','k','Linewidth',2,'Parent',h3); hold on;
%print('-dpng',[ppath 'Hist_Fore_',harv,'_tot_fracs_FD_logZlDet_temp_diff1951.png'])

%%
figure
line(time,plank,'color',cmap_ppt(5,:),'Linewidth',2); hold on;
line(time,oxy,'color',cmap_ppt(1,:),'Linewidth',2); hold on;
line(time,temp2,'color',cm{6},'Linewidth',2); hold on;
line(time,temp1,'color','k','Linewidth',2); hold on;
line(time,temp3,'color',cmap_ppt(3,:),'Linewidth',2); hold on;

%%
figure(5)
% line(time,plank,'color',cmap_ppt(5,:),'Linewidth',2); hold on;
% line(time,eplank,'color',cmap_ppt(5,:)); hold on;
line(time,plank,'color',cm{16},'Linewidth',2); hold on;
line(time,eplank,'color',cm{16}); hold on;
xlim([time(1) time(end)])
set(gca,'XTickLabels',[],'YTickLabels',[])
xlabel('Time')
title('Plankton')
print('-dpng',[ppath 'Ensemble_plankton.png'])

figure(6)
% line(time,oxy,'color',cmap_ppt(1,:),'Linewidth',2); hold on;
% line(time,eoxy,'color',cmap_ppt(1,:)); hold on;
line(time,oxy,'color',cm{18},'Linewidth',2); hold on;
line(time,eoxy,'color',cm{18}); hold on;
xlim([time(1) time(end)])
set(gca,'XTickLabels',[],'YTickLabels',[])
xlabel('Time')
title('Oxygen')
%cm{17}
%cm{18}
print('-dpng',[ppath 'Ensemble_oxygen.png'])

figure(7)
line(time,temp1,'color','k','Linewidth',2); hold on;
line(time,etemp1,'color','k'); hold on;
xlim([time(1) time(end)])
set(gca,'XTickLabels',[],'YTickLabels',[])
xlabel('Time')
title('Temperature')
print('-dpng',[ppath 'Ensemble_temperature.png'])

%%
figure(8)
line(time,temp2,'color',cm{6},'Linewidth',2); hold on;
line(time,etemp2,'color',cm{6}); hold on;

figure(9)
line(time,temp3,'color',cmap_ppt(3,:),'Linewidth',2); hold on;
line(time,etemp3,'color',cmap_ppt(3,:)); hold on;


