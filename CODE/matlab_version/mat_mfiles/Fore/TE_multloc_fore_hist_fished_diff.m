% Changes in TE from Hist 1951-2000 to Fore 2051-2100

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = [dp cfile '/'];
ppath = [pp cfile '/'];

% Hist grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([gpath 'grid_csv.csv']);
ID = grid(:,1);

geolon_t=double(geolon_t);
geolat_t=double(geolat_t);
[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% FEISTY Hist 1951-2000
load([fpath 'TEeffDet_Historic_All_fish03_' cfile '.mat'],'TEeff_L',...
    'TEeff_HTLd','TEeff_LTLd');

hLe = TEeff_L;
hHTLe = TEeff_HTLd;
hLTLe = TEeff_LTLd;

clear TEeff_L TEeff_HTLd TEeff_LTLd

%% FEISTY Fore 2051-2100
load([fpath 'TEeffDet_Forecast_All_fish03_' cfile '.mat'],'TEeff_L',...
    'TEeff_HTLd','TEeff_LTLd');

fLe = TEeff_L;
fHTLe = TEeff_HTLd;
fLTLe = TEeff_LTLd;

clear TEeff_L TEeff_HTLd TEeff_LTLd

%% calc TEs from effective TEs
hLTL = real(hLTLe.^(1/(4/3)));
hATL = real(hLe.^(1/4));         
hHTL = real(hHTLe.^(1/3));

fLTL = real(fLTLe.^(1/(4/3)));
fATL = real(fLe.^(1/4));         
fHTL = real(fHTLe.^(1/3));

%% Percent diffs in TEeff
pdiffL = (fLe-hLe) ./ hLe;
pdiffHTL = (fHTLe-hHTLe) ./ hHTLe;
pdiffLTL = (fLTLe-hLTLe) ./ hLTLe;

%Pure diff in TE
diffL = (fATL-hATL) ;
diffHTL = (fHTL-hHTL) ;
diffLTL = (fLTL-hLTL) ;

%% Pdiff in TEeff
figure(1)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pdiffLTL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Percent difference TEeff LTL')

%HTL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pdiffHTL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('Percent difference TEeff HTL')

%L
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,pdiffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Percent difference TEeff ATL')
%stamp(cfile)
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_TEeff_pdiff_subplot.png'])

%% Diff in TE
figure(2)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffLTL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.1 0.1]);
set(gcf,'renderer','painters')
title('Difference TE LTL')

%HTL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffHTL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.1 0.1]);
set(gcf,'renderer','painters')
title('Difference TE HTL')

%L
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.1 0.1]);
colorbar('Position',[0.2 0.495 0.6 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Difference TE ATL')
%stamp(cfile)
print('-dpng',[ppath 'Hist_Fore_',harv,'_global_TE_diff_subplot.png'])

%%
% Individual Hist vs Fore
figure(3)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hLTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.25]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast TE LTL');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,fLTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.25]);
colorbar
set(gcf,'renderer','painters')
title('Forecast TE LTL');
print('-dpng',[ppath 'Hist_Fore_' harv '_global_TE_LTL.png'])

%%
figure(4)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hHTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.4]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast TE HTL');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,fHTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.4]);
colorbar
set(gcf,'renderer','painters')
title('Forecast TE HTL');
print('-dpng',[ppath 'Hist_Fore_' harv '_global_TE_HTL.png'])

%%
figure(5)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hATL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
colorbar
set(gcf,'renderer','painters')
title('Hindcast TE ATL');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,fATL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
colorbar
set(gcf,'renderer','painters')
title('Forecast TE ATL');
print('-dpng',[ppath 'Hist_Fore_' harv '_global_TE_ATL.png'])

%% Subplot of historic TEs
figure(6)
subplot('Position',[0 0.53 0.5 0.5])
%LTL
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hLTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
set(gcf,'renderer','painters')
title('TE LTL')

%HTL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hHTL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
set(gcf,'renderer','painters')
title('TE HTL')

%L
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,hATL)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.05 0.3]);
colorbar('Position',[0.2 0.495 0.6 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('TE ATL')
%stamp(cfile)
print('-dpng',[ppath 'Hist_',harv,'_global_TE_subplot.png'])

