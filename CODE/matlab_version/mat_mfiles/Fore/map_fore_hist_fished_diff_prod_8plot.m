% Visualize difference between
% ESM2M Hindcast of 1951-2000 
% and Forecast of 2051-2100
% PRODUCTION

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
dpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
load([dpath 'Prod_pdiffs_hist_fore_',harv,'_' cfile '.mat'],...
    'pdiffF','pdiffP','pdiffD','pdiffB','pdiffAll',...
    'pdiffN','pdiffDet','pdiffMZ','pdiffLZ','pdiffZ','diffZD')

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% figure info
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - npp
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffN)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Net primary','HorizontalAlignment','center')
text(-2.5,1.75,'A')

% - mesozoo
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Mesozoo','HorizontalAlignment','center')
text(-2.5,1.75,'B')

% - det
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Detritus','HorizontalAlignment','center')
text(-2.5,1.75,'C')

% - benthos
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.5,1.75,'D')

% - forage
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.5,1.75,'E')

% - LP
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'F')

%F - dem
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'G')

%- All
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'H')

print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_prod_8plot.png'])

%% Mollweid proj

f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - npp
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffN)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Net primary','HorizontalAlignment','center')
text(-2.5,1.75,'A')

%B - mesozoo
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffZ)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Mesozoo','HorizontalAlignment','center')
text(-2.5,1.75,'B')

%C - det
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffDet)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Detritus','HorizontalAlignment','center')
text(-2.5,1.75,'C')

%D - benthos
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffB)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.5,1.75,'D')

%E - forage
subplot('Position',[0.475 0.75 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.5,1.75,'E')

%F - LP
subplot('Position',[0.475 0.5 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.5,1.75,'F')

%G - dem
subplot('Position',[0.475 0.25 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.5,1.75,'G')

%H - All
subplot('Position',[0.475 0.0 0.44 0.25])
axesm ('Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,100*pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
text(-2.5,1.75,'H')

print('-dpng',[pp 'Hist_Fore_',harv,'_global_pdiff_prod_8plot_mollw.png'])