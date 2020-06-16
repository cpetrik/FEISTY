% Map global benthos with Climatol and Wei et al

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
datap = '/Volumes/FEISTY/NC/Matlab_new_size/';
cp = '/Volumes/FEISTY/CSV/Matlab_new_size/';

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cdir='/Volumes/FEISTY/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% colors
cmBP=cbrewer('seq','BuPu',50,'PCHIP');

%% Model bent & dem
[ni,nj]=size(lon);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_2B75_BE08_noCC_RE00100';
dpath = [datap cfile '/'];
fpath = [figp cfile '/'];
load([dpath 'Means_Climatol_All_fish03_' cfile '.mat'],...
    'md_mean','ld_mean','sb_mean','mb_mean');
Zmd8 = NaN*ones(ni,nj);
Zld8 = NaN*ones(ni,nj);
Zsb8 = NaN*ones(ni,nj);
Zmb8 = NaN*ones(ni,nj);
Zmd8(ID) = md_mean;
Zld8(ID) = ld_mean;
Zsb8(ID) = sb_mean;
Zmb8(ID) = mb_mean;
Zb8 = Zsb8 + Zmb8;
Zd8 = Zmd8 + Zld8;
clear md_mean ld_mean b_mean

%% Sizes
%FEISTY small bent 4.6?36.8 (13)                4.6?37 mm
%FEISTY medium bent/dem 36.8?292.4 (104)        37?292 mm
%FEISTY large dem 292.4?2320.8 (824)            292?2321 mm                         
%Wei meio infauna 20-74 um =                    0.020-0.074 mm
%Wei macro infauna 250-520 um =                 0.25-0.52 mm
%Wei mega = epibenthic inverts + fish >1 cm =   >10 mm
%Wei invert >1 cm =                             >10 mm
%Wei fish >1 cm =                               >10 mm

%% Wei benthic biomass
seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei2010_Global_seafloor_biomass.csv',1,0);
Wcol = {'latitude','longitude','depth','bact.biom.mean','meio.biom.mean',...
    'macro.biom.mean','mega.biom.mean','inv.biom.mean','fis.biom.mean'};
Wcol = Wcol';

% all mean biomasses in log10 mg C/m2
meio = seafl(:,5);
macro = seafl(:,6);
mega = seafl(:,7);
invert = seafl(:,8);
fish = seafl(:,9);
% convert to g WW/m2
meio = 10.^(meio) * 1e-3 * 9.0;
macro = 10.^(macro) * 1e-3 * 9.0;
mega = 10.^(mega) * 1e-3 * 9.0;
invert = 10.^(invert) * 1e-3 * 9.0;
fish = 10.^(fish) * 1e-3 * 9.0;

%% put on same grid as POEM output
land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

% SHIFT
wlat = seafl(:,1);
wlon = seafl(:,2);

%ESM lon 0.5 = -100
%100.5:179.5 needs to be 0.5:80.5
test=wlon;
id=(wlon>=100);
test(id==1)=wlon(id==1)-100;
%-179:99 needs to be 81.5:359.5
test(id==0)=wlon(id==0)+260;
tlon1 = wlon(id==1)-100;
tlon2 = wlon(id==0)+260;
swlon = test;

[geolon_w,geolat_w] = meshgrid(-179.5:179.5,-89.5:89.5);

Zmi = griddata(wlat,wlon,meio,geolat_w,geolon_w);
Zma = griddata(wlat,wlon,macro,geolat_w,geolon_w);
Zme = griddata(wlat,wlon,mega,geolat_w,geolon_w);
Zi = griddata(wlat,wlon,invert,geolat_w,geolon_w);
Zf = griddata(wlat,wlon,fish,geolat_w,geolon_w);

fauna = Zmi + Zma + Zme;
ifauna = Zmi + Zma;

%% Interpolate to same grid
%lat        [-89.5 89.5]
%lon        [0.5 359.5]

% Need to fix longitude
test = lon-360;
id=find(test<-180);
test(id)=test(id)+360;
lon = test;

lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

clat = lat;
clon = lon;

Zsb = griddata(clat,clon,Zsb8,glat,glon);
Zmb = griddata(clat,clon,Zmb8,glat,glon);
Zb = griddata(clat,clon,Zb8,glat,glon);
Zd = griddata(clat,clon,Zd8,glat,glon);

% plot info
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%% Maps
figure(1)
% meio
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zmi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
colorbar
set(gcf,'renderer','painters')
title('Wei log_1_0 meio-infauna (g m^-^2)')

subplot('Position',[0.5 0.51 0.5 0.5])
% macro
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zma))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
colorbar
set(gcf,'renderer','painters')
title('Wei log_1_0 macro-infauna (g m^-^2)')

subplot('Position',[0 0 0.5 0.5])
% mega
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zme))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
colorbar
set(gcf,'renderer','painters')
title('Wei log_1_0 mega-epifauna (g m^-^2)')
print('-dpng',[fpath 'Global_invert_fauna_Wei.png'])

%%
figure(2)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zma))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%axesmui
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean macro-infauna (g m^-^2)')

% FEISTY
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 all benthos (g m^-^2)')

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zsb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 small benthos (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zmb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 medium benthos (g m^-^2)')
print('-dpng',[fpath 'Clim_All_fish03_global_macrofauna_Wei.png'])

%%
figure(3)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zme))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%axesmui
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. mean megafauna (g m^-^2)')

% FEISTY
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY all benthos (g m^-^2)')

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zsb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY small benthos (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zmb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY medium benthos (g m^-^2)')
print('-dpng',[fpath 'Clim_All_fish03_global_megafauna_Wei.png'])

%%
figure(4)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zma+Zi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%axesmui
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. mean macro and mega invertebrates')

% FEISTY
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY all benthos (g m^-^2)')

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zsb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY small benthos (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zmb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY medium benthos (g m^-^2)')
print('-dpng',[fpath 'Clim_All_fish03_global_macro_inverts_Wei.png'])

%%
figure(4)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1);%,'origin',[0 -100 0] ,'FLonLimit',[-180 180])
surfm(geolat_w,geolon_w,(Zma+Zi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%axesmui
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. mean macro and mega invertebrates')

% FEISTY
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1);%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY all benthos (g m^-^2)')

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zsb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY small benthos (g m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(Zmb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e2]);
set(gcf,'renderer','painters')
title('FEISTY medium benthos (g m^-^2)')
print('-dpng',[fpath 'Clim_All_fish03_global_macro_inverts_Wei.png'])

%% Put on same axes and take difference
% comp small w/ all infauna, ifauna
% med w/ all invert epifauna (no fish), Zi
% all w/ all invert (no fish), fauna

diffS = (Zsb - ifauna);
diffM = (Zmb - ifauna);
diffB = (Zb - fauna);

figure(5)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1);%,'origin',[0 -100 0] ,'FLonLimit',[-180 180])
surfm(geolat_w,geolon_w,(diffS))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e3]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Diff Wei infauna from S benthos (g m^-^2)')

% FEISTY
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1);%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,(diffM))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e3]);
set(gcf,'renderer','painters')
title('Diff Wei epi inverts from M benthos (g m^-^2)')

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_w,geolon_w,(diffB))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gca,'ColorScale','log')
caxis([1e-2 1e3]);
set(gcf,'renderer','painters')
title('Diff Wei all inverts from All benthos (g m^-^2)')

% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_w,geolon_w,(Zmb))
% colormap(cmBP)
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1.5 2]);
% set(gcf,'renderer','painters')
% title('FEISTY medium benthos (g m^-^2)')
print('-dpng',[fpath 'Clim_All_fish03_global_inverts_Wei_diff.png'])



