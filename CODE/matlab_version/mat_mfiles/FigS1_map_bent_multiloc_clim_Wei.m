% Map global benthos with Climatol and Wei et al

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
datap = '/Volumes/GFDL/NC/Matlab_new_size/';
cp = '/Volumes/GFDL/CSV/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% colors
cmBP=cbrewer('seq','BuPu',50);

%% Model bent & dem
[ni,nj]=size(lon);

cfile = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_',...
    'BE08_noCC_RE00100'];
dpath = [datap cfile '/'];
fpath = [figp cfile '/'];
load([dpath 'Means_Climatol_All_fish03_' cfile '.mat'],...
    'md_mean','ld_mean','b_mean');
Zmd8=NaN*ones(ni,nj);
Zld8=NaN*ones(ni,nj);
Zb8 =NaN*ones(ni,nj);
Zmd8(ID)=md_mean;
Zld8(ID)=ld_mean;
Zb8(ID) =b_mean;
Zd8 = Zmd8+Zld8;
clear md_mean ld_mean b_mean

cfile2 = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_',...
    'BE03_noCC_RE00100'];
dpath2 = [datap cfile2 '/'];
load([dpath2 'Climatol_All_fish03.mat'],...
    'md_mean','ld_mean','b_mean');
Zmd3=NaN*ones(ni,nj);
Zld3=NaN*ones(ni,nj);
Zb3 =NaN*ones(ni,nj);
Zmd3(ID)=md_mean;
Zld3(ID)=ld_mean;
Zb3(ID)=b_mean;
Zd3 = Zmd3+Zld3;
clear md_mean ld_mean b_mean

cfile3 = ['Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_',...
    'BE13_noCC_RE00100'];
dpath3 = [datap cfile3 '/'];
load([dpath3 'Climatol_All_fish03.mat'],...
    'md_mean','ld_mean','b_mean');

Zmd13=NaN*ones(ni,nj);
Zld13=NaN*ones(ni,nj);
Zb13=NaN*ones(ni,nj);
Zmd13(ID)=md_mean;
Zld13(ID)=ld_mean;
Zb13(ID)=b_mean;
Zd13 = Zmd13+Zld13;
clear md_mean ld_mean b_mean

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
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

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


%% Maps
figure(1)
% meio
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zmi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
colorbar
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean meiofauna (g m^-^2)')

figure(2)
% macro
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zma))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
colorbar
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean macrofauna (g m^-^2)')

figure(3)
% mega
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zme))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
colorbar
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean megafauna (g m^-^2)')


%% Inv
figure(4)
% Wei
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_w,geolon_w,log10(Zi))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean invertebrates (g m^-^2)')

% FEISTY
subplot(2,2,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb3))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.025')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb8))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.075')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb13))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.125')
%print('-dpng',[fpath 'Clim_All_fish03_global_benthos_Wei.png'])

%%
figure(5)
% Wei
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(Zi+Zf))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-1.5 2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean benthos (g m^-^2)')

% FEISTY
subplot(2,2,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb3))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.025')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb8))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.075')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb13))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.125')
%print('-dpng',[fpath 'Clim_All_fish03_global_benthos_Wei.png'])

%%
fauna = Zmi + Zme + Zma;
efauna = Zme + Zma;
figure(6)
% Wei
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(fauna))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-1.5 2]);
colorbar('Position',[0.05 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean benthos (g m^-^2)')

% FEISTY
subplot(2,2,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb3))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.025')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb8))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.075')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb13))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.125')
%print('-dpng',[fpath 'Clim_All_fish03_global_fauna_Wei.png'])

%%
figure(7)
% Wei
subplot(2,2,1)
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
subplot(2,2,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb3))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.025')

subplot(2,2,3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb8))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.075')

subplot(2,2,4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb13))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.125')
%print('-dpng',[fpath 'Clim_All_fish03_global_epifauna_Wei.png'])

%%
figure(8)
% Wei
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0] ,'FLonLimit',[-180 180])% thos doesn't work for some reason
surfm(geolat_w,geolon_w,log10(Zma))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
axesmui
caxis([-1.5 2]);
colorbar('Position',[0.25 0.525 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Wei et al. log_1_0 mean macrofauna (g m^-^2)')

% FEISTY
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb3))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.025')

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb8))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.075')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb13))
colormap(cmBP)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2]);
set(gcf,'renderer','painters')
title('FEISTY log_1_0 mean benthos (g m^-^2), \beta=0.125')
print('-dpng',[fpath 'Clim_All_fish03_global_macrofauna_Wei_multBE.png'])

