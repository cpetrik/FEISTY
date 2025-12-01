% Compare pelagic and demersal biomass or catch
% to observed effort
% spatially

clear 
close all

%% Effort data 
epath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/movement/FishingDataColleen/';
load([epath 'Global_fishing_effort_1deg_pel_dem_mean2015-2024.mat']);

%grid lat & lon
[elat,elon] = meshgrid(lat,lon);

%%
figure(1)
pcolor(meanPeffort); shading flat
colorbar
clim([0 10])

figure(2)
pcolor(elat); shading flat
colorbar
title('lat')

figure(3)
pcolor(elon); shading flat
colorbar
title('lon')

%% FEISTY
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%%
figure(4)
pcolor(geolat_t); shading flat
colorbar
title('Glat')

figure(5)
pcolor(geolon_t); shading flat
colorbar
title('Glon')

%%
exper = 'CORE_Hindcast1988_no_move_';
load([fpath 'Means_' exper cfile '.mat']);

% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(GRD.ID)=sf_mean20;
Zsp(GRD.ID)=sp_mean20;
Zsd(GRD.ID)=sd_mean20;
Zmf(GRD.ID)=mf_mean20;
Zmp(GRD.ID)=mp_mean20;
Zmd(GRD.ID)=md_mean20;
Zlp(GRD.ID)=lp_mean20;
Zld(GRD.ID)=ld_mean20;
Zb(GRD.ID)=b_mean20;

sAll = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
sAllF = Zsf+Zmf;
sAllP = Zsp+Zmp+Zlp;
sAllD = Zsd+Zmd+Zld;

%% Interp to same grid
%elat       [-89.7500 89.2500]
%geolat_t   [-81.5 89.4879]
%elon       [-179.7500 179.2500]
%geolon_t   [-279.9803 79.9803]

%% Need to fix GFDL longitude
test2=geolon_t;
id=find(test2<-180);
test2(id)=test2(id)+360;
geolon = test2;

geolat = geolat_t;
geolon = geolon;

%%
figure(6)
pcolor(geolat); shading flat
colorbar
title('Glat2')

figure(7)
pcolor(geolon); shading flat
colorbar
title('Glon2')

%%
lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

%%
cAll = griddata(geolat,geolon,sAll,glat,glon);
cP = griddata(geolat,geolon,sAllP,glat,glon);
cD = griddata(geolat,geolon,sAllD,glat,glon);
cF = griddata(geolat,geolon,sAllF,glat,glon);

eP = griddata(elat,elon,meanPeffort,glat,glon);
eD = griddata(elat,elon,meanDeffort,glat,glon);

%% Pelagic
figure(3)
subplot('Position',[0 0.53 0.5 0.5])
%No move
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,cF)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('Forage')

%Move
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,cP)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('Lg Pel')

%Diff
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,eP)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-0.01 0.01]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Pel effort')
%print('-dpng',[ppath 'COREstationary_moveprey_global_BENT.png'])

