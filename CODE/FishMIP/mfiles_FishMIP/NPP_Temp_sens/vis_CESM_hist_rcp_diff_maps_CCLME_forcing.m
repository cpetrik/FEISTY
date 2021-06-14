% Visualize output from CESM
% CCLME maps of 2090-2100 vs 1990-2000

clear all
close all

%% Map data
cpath = '/Volumes/FEISTY/Fish-MIP/CMIP5/CESM/';
load([cpath 'gridspec_cesm.mat']);
load([cpath 'Data_grid_cesm.mat']);
[ni,nj]=size(LON);
% plotminlat=-90; 
% plotmaxlat=90;
% plotminlon=-280;
% plotmaxlon=80;
plotminlat=25; 
plotmaxlat=55;
plotminlon=-134;
plotmaxlon=-110;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines

%% CESM Hist data
ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
load([cpath 'Hist/cesm_hist_tp_100_monthly_185001-200512.mat'])
load([cpath 'Hist/cesm_hist_szoo_100_monthly_185001-200512.mat'])
load([cpath 'Hist/cesm_hist_lzoo_100_monthly_185001-200512.mat'])
time = (1850+(15/365)):(1/12):(2006);

%% 
%Comparing 2090-2100 to 1990-2000
Hyr = time;
yrP=find(Hyr>1990 & Hyr<=2000);  

tp100_P = double(mean(tp_100(:,:,yrP),3));
sz100_P = double(mean(nsmz_100(:,:,yrP),3));
lz100_P = double(mean(nlgz_100(:,:,yrP),3));

clear tp_100 nlgz_100 nsmz_100

%% CESM Fore data
load([cpath 'RCP85/cesm_rcp85_tp_100_monthly_200601-210012.mat'])
load([cpath 'RCP85/cesm_rcp85_szoo_100_monthly_200601-210012.mat'])
load([cpath 'RCP85/cesm_rcp85_lzoo_100_monthly_200601-210012.mat'])

time = (2006+(15/365)):(1/12):(2101);

% Comparing 2090-2100 to 1990-2000
Fyr = time;
yrF=find(Fyr>2090 & Fyr<=2100); 

tp100_F = double(mean(tp_100(:,:,yrF),3));
sz100_F = double(mean(nsmz_100(:,:,yrF),3));
lz100_F = double(mean(nlgz_100(:,:,yrF),3));

clear tp_100 nlgz_100 nsmz_100

%% Diff maps of all fish
PFracLM = lz100_P ./ (lz100_P + sz100_P);
FFracLM = lz100_F ./ (lz100_F + sz100_F);

diffTP = (tp100_F - tp100_P);
diffLM = (FFracLM - PFracLM);

pdiffS = (sz100_F - sz100_P) ./ sz100_P;
pdiffL = (lz100_F - lz100_P) ./ lz100_P;

%% All 4 on subplots 
figure(1)
% TP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffTP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
%colorbar('Position',[0.475 0.25 0.03 0.5],'orientation','vertical')
colorbar
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Upper ocean temperature','HorizontalAlignment','center')

% O2
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,PFracLM)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf O_2','HorizontalAlignment','center')

% pH
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,FFracLM)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf pH','HorizontalAlignment','center')

% plankton size
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffLM)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.2 0.2]);
colorbar
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Fraction large phyto','HorizontalAlignment','center')
%print('-dpng',[ppath 'Hist_RCP_CCLME_bgc_subplot_diffs.png'])

%% All 4 on subplots 
figure(2)
% SmZ
subplot('Position',[0.05 0.51 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffS)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-30 30]);
colorbar('Position',[0.325 0.55 0.03 0.4],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Small plankton','HorizontalAlignment','center')

% LgZ
subplot('Position',[0.375 0.51 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,100*pdiffL)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-30 30]);
colorbar('Position',[0.65 0.55 0.03 0.4],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Large plankton','HorizontalAlignment','center')
text(-0.003,0.95,'RCP 8.5 2090-2100 vs ','HorizontalAlign','left')
text(-0.003,0.90,'Historic 1990-2000','HorizontalAlign','left')
text(-0.003,0.81,'% change in biomass','HorizontalAlign','left')

% plankton size
subplot('Position',[0.05 0 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffLM)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.2 0.2]);
colorbar('Position',[0.325 0.04 0.03 0.4],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Fraction large plankton','HorizontalAlignment','center')

% TP
subplot('Position',[0.375 0 0.3 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,diffTP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
colorbar('Position',[0.65 0.04 0.03 0.4],'orientation','vertical')
set(gcf,'renderer','painters')
text(-0.30,0.935,'\bf Upper ocean temperature','HorizontalAlignment','center')
text(-0.003,0.91,'Absolute change','HorizontalAlign','left')
print('-dpng',[ppath 'Hist_RCP_CCLME_bgc_subplot_diffs.png'])

