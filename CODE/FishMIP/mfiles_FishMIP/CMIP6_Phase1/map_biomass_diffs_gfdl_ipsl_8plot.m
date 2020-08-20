% Visualize difference between Pre, Hist, SSP126, SSP585
% GFDL & IPSL models

clear all
close all

%% Fish data
cfile = 'Dc_Lam579_enc70-b200_m440-b175-k086_c20-b250_D080_A050_nmort1_BE10_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% FEISTY Output
%gfdl
gpath=['/Volumes/FEISTY/NC/FishMIP/GFDL_CMIP6/' cfile '/'];
% load([gpath 'Means_Pre_' cfile '.mat'],...
%     'GPreAllF','GPreAllP','GPreAllD','GPreAllM','GPreAllL','GPreAll')
load([gpath 'Means_Hist_' cfile '.mat'],...
    'GHistAllF','GHistAllP','GHistAllD','GHistAllM','GHistAllL','GHistAll')
load([gpath 'Means_SSP126_2090-2100_' cfile '.mat'],...
    'GS126AllF','GS126AllP','GS126AllD','GS126AllM','GS126AllL','GS126All')
load([gpath 'Means_SSP585_2090-2100_' cfile '.mat'],...
    'GS585AllF','GS585AllP','GS585AllD','GS585AllM','GS585AllL','GS585All');

%ipsl
ipath=['/Volumes/FEISTY/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
% load([ipath 'Means_Pre_' cfile '.mat'],...
%     'IPreAllF','IPreAllP','IPreAllD','IPreAllM','IPreAllL','IPreAll')
load([ipath 'Means_Hist_' cfile '.mat'],...
    'IHistAllF','IHistAllP','IHistAllD','IHistAllM','IHistAllL','IHistAll')
load([ipath 'Means_SSP126_2090-2100_' cfile '.mat'],...
    'IS126AllF','IS126AllP','IS126AllD','IS126AllM','IS126AllL','IS126All')
load([ipath 'Means_SSP585_2090-2100_' cfile '.mat'],...
    'IS585AllF','IS585AllP','IS585AllD','IS585AllM','IS585AllL','IS585All');

%% Grid info
cpath = '/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/';
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/gridspec_gfdl_cmip6.mat');
load([cpath 'Data_grid_gfdl.mat']);

[ni,nj]=size(LON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% Diffs

GdiffH126_F = (GS126AllF - GHistAllF) ./ GHistAllF;
GdiffH585_F = (GS585AllF - GHistAllF) ./ GHistAllF;
GdiffH126_P = (GS126AllP - GHistAllP) ./ GHistAllP;
GdiffH585_P = (GS585AllP - GHistAllP) ./ GHistAllP;
GdiffH126_D = (GS126AllD - GHistAllD) ./ GHistAllD;
GdiffH585_D = (GS585AllD - GHistAllD) ./ GHistAllD;
GdiffH126_A = (GS126All  - GHistAll)  ./ GHistAll;
GdiffH585_A = (GS585All  - GHistAll)  ./ GHistAll;

IdiffH126_F = (IS126AllF - IHistAllF) ./ IHistAllF;
IdiffH585_F = (IS585AllF - IHistAllF) ./ IHistAllF;
IdiffH126_P = (IS126AllP - IHistAllP) ./ IHistAllP;
IdiffH585_P = (IS585AllP - IHistAllP) ./ IHistAllP;
IdiffH126_D = (IS126AllD - IHistAllD) ./ IHistAllD;
IdiffH585_D = (IS585AllD - IHistAllD) ./ IHistAllD;
IdiffH126_A = (IS126All  - IHistAll)  ./ IHistAll;
IdiffH585_A = (IS585All  - IHistAll)  ./ IHistAll;

%% figure info
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - GFDL SSP126-Hist
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH126_F)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'GFDL')

%B - GFDL SSP585-Hist
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH585_F)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %C - GFDL SSP126-Pre
% subplot('Position',[0.025 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP126_F)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %D - GFDL SSP585-Pre
% subplot('Position',[0.025 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP585_F)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')


%E - IPSL SSP126-Hist
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH126_F)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'IPSL')

%F - IPSL SSP585-Hist
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH585_F)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %G - IPSL SSP126-Pre
% subplot('Position',[0.475 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP126_F)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %H - IPSL SSP585-Pre
% subplot('Position',[0.475 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP585_F)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')

print('-dpng',[pp '_Pre_Hist_SSPs_global_pdiff_biom_8plot_F.png'])

%% Large pel
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f2.Units = 'inches';

%A - GFDL SSP126-Hist
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH126_P)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'GFDL')

%B - GFDL SSP585-Hist
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH585_P)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %C - GFDL SSP126-Pre
% subplot('Position',[0.025 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP126_P)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %D - GFDL SSP585-Pre
% subplot('Position',[0.025 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP585_P)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')


%E - IPSL SSP126-Hist
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH126_P)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'IPSL')

%F - IPSL SSP585-Hist
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH585_P)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %G - IPSL SSP126-Pre
% subplot('Position',[0.475 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP126_P)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %H - IPSL SSP585-Pre
% subplot('Position',[0.475 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP585_P)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')

print('-dpng',[pp '_Pre_Hist_SSPs_global_pdiff_biom_8plot_P.png'])

%% Demersal
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f3.Units = 'inches';

%A - GFDL SSP126-Hist
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH126_D)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'GFDL')

%B - GFDL SSP585-Hist
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH585_D)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %C - GFDL SSP126-Pre
% subplot('Position',[0.025 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP126_D)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %D - GFDL SSP585-Pre
% subplot('Position',[0.025 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP585_D)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')


%E - IPSL SSP126-Hist
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH126_D)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'IPSL')

%F - IPSL SSP585-Hist
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH585_D)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %G - IPSL SSP126-Pre
% subplot('Position',[0.475 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP126_D)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %H - IPSL SSP585-Pre
% subplot('Position',[0.475 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP585_D)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')

print('-dpng',[pp '_Pre_Hist_SSPs_global_pdiff_biom_8plot_D.png'])

%% All fish
f4 = figure('Units','inches','Position',[1 3 6.5 8]);
%f4.Units = 'inches';

%A - GFDL SSP126-Hist
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH126_A)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'GFDL')

%B - GFDL SSP585-Hist
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*GdiffH585_A)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %C - GFDL SSP126-Pre
% subplot('Position',[0.025 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP126_A)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %D - GFDL SSP585-Pre
% subplot('Position',[0.025 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*GdiffP585_A)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')


%E - IPSL SSP126-Hist
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH126_A)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'SSP 126 - Hist','HorizontalAlignment','center')
text(-2.5,1.75,'IPSL')

%F - IPSL SSP585-Hist
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,100*IdiffH585_A)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SSP 585 - Hist','HorizontalAlignment','center')

% %G - IPSL SSP126-Pre
% subplot('Position',[0.475 0.25 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP126_A)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 126 - Pre','HorizontalAlignment','center')
% 
% %H - IPSL SSP585-Pre
% subplot('Position',[0.475 0.0 0.4 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(LAT,LON,100*IdiffP585_A)
% cmocean('balance')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-100 100]);
% set(gcf,'renderer','painters')
% text(0,1.75,'SSP 585 - Pre','HorizontalAlignment','center')

print('-dpng',[pp '_Pre_Hist_SSPs_global_pdiff_biom_8plot_All.png'])


