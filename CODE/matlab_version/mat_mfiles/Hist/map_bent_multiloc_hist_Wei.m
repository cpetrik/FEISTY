% Map global benthos with Climatol and Wei et al

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
datap = '/Volumes/GFDL/NC/Matlab_new_size/';
cp = '/Volumes/GFDL/CSV/Matlab_new_size/';

% colors
cmBP=cbrewer('seq','BuPu',50);

cfile = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_',...
    'BE08_noCC_RE00100'];
harv = 'All_fish03';
dpath = [datap cfile '/'];
fpath = [figp cfile '/'];

load([dpath 'Means_Historic_',harv,'_' cfile '.mat'],...
    'md_mean5','ld_mean5','b_mean5');

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);


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
% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

[ni,nj]=size(geolon_t);
ID = grid(:,1);
land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%% SHIFTED WRONG OR SOMETHING
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

fauna = Zmi + Zme + Zma;
efauna = Zme + Zma;

%% Model bent & dem
Zmd=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zmd(ID)=md_mean5;
Zld(ID)=ld_mean5;
Zb(ID)=b_mean5;

Zd = Zmd+Zld;

%% Maps
figure(10)
% meio
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zmi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean meiofauna (g m^-^2)')

figure(11)
% macro
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zma))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean macrofauna (g m^-^2)')

figure(12)
% mega
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zme))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean megafauna (g m^-^2)')


%% Inv
figure(1)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.01 0.525 0.48 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean invertebrates (g m^-^2)')

% FEISTY
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2)')
print('-dpng',[fpath 'Hist9095_All_fish03_global_inverts_Wei.png'])

%% Fish
figure(2)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zf))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.01 0.525 0.48 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean fishes (g m^-^2)')

% FEISTY
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zd))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean demersals (g m^-^2)')
print('-dpng',[fpath 'Hist9095_All_fish03_global_Demersal_Wei.png'])

%%
figure(3)
% Wei
subplot('Position',[0.05 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(Zi+Zf))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-2 1]);
colorbar('Position',[0.05 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean inverts & fish (g m^-^2)')

% FEISTY
subplot('Position',[0.05 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2)')
print('-dpng',[fpath 'Hist9095_All_fish03_global_benthos_Wei.png'])

%%
figure(4)
% Wei
subplot('Position',[0.05 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(fauna))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-1 2]);
colorbar('Position',[0.05 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean infauna & epifauna (g m^-^2)')

% FEISTY
subplot('Position',[0.05 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2)')
print('-dpng',[fpath 'Hist9095_All_fish03_global_fauna_Wei.png'])

%%
figure(5)
% Wei
subplot('Position',[0.05 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(efauna))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-1.5 2]);
colorbar('Position',[0.05 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean epifauna (g m^-^2)')

% FEISTY
subplot('Position',[0.05 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2)')
print('-dpng',[fpath 'Hist9095_All_fish03_global_epifauna_Wei.png'])

%%
figure(6)
% Wei
subplot('Position',[0.05 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(Zma))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-2 2]);
colorbar('Position',[0.05 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean macrofauna (g m^-^2)')

% FEISTY
subplot('Position',[0.05 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2)')
print('-dpng',[fpath 'Hist9095_All_fish03_global_macrofauna_Wei.png'])

