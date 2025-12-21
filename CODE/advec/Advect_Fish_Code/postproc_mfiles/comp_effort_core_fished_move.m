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
exper = 'CORE_Hindcast1988_no_move_All_fish03_yield_';
load([fpath 'Annual_Means_' exper cfile '.mat']);

CP = Cmp + Clp;
CD = Cmd + Cld;
CPel = Cmf + CP;
CAll = CPel + CD;

% Plots in space
ZF = mean(Cmf,3,'omitnan');
ZP = mean(CP,3,'omitnan');
ZD = mean(CD,3,'omitnan');
ZPel = mean(CPel,3,'omitnan');
ZAll = mean(CAll,3,'omitnan');

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
cAll = griddata(geolat,geolon,ZAll,glat,glon);
cPel = griddata(geolat,geolon,ZPel,glat,glon);
cP = griddata(geolat,geolon,ZP,glat,glon);
cD = griddata(geolat,geolon,ZD,glat,glon);
cF = griddata(geolat,geolon,ZF,glat,glon);

eP = griddata(elat,elon,meanPeffort,glat,glon);
eD = griddata(elat,elon,meanDeffort,glat,glon);

cP(cP<=0) = nan;

%% Pelagic - check things aligned
figure(3)
subplot('Position',[0 0.53 0.49 0.49])
%F
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,log10(cF))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([6 12]);
colorbar
set(gcf,'renderer','painters')
title('Forage')

%LP
subplot('Position',[0.5 0.53 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,log10(cP))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([6 12]);
colorbar
set(gcf,'renderer','painters')
title('Lg Pel')

%Effort
subplot('Position',[0.25 0.0 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,log10(eP))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Pel effort')
print('-dpng',[ppath exper 'effort_global_Pelagics.png'])


%% Demersal - check things aligned
figure(4)
subplot('Position',[0 0.53 0.49 0.49])
%F
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,log10(cD))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([6 12]);
colorbar
set(gcf,'renderer','painters')
title('Demersal')

%LP
% subplot('Position',[0.5 0.53 0.49 0.49])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(glat,glon,log10(cP))
% cmocean('dense')
% load coastlines;                     
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([6 12]);
% colorbar
% set(gcf,'renderer','painters')
% title('Lg Pel')

%Effort
subplot('Position',[0.25 0.0 0.49 0.49])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,log10(eD))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('Dem effort')
print('-dpng',[ppath exper 'effort_global_Demersals.png'])

%% Scatter plots
figure(5)
subplot(2,2,1)
plot(log10(eP(:)+1),log10(cPel(:)+1),'.k')
title('All Pel')
ylim([6 12])

subplot(2,2,2)
plot(log10(eP(:)+1),log10(cF(:)+1),'.r')
title('F')
ylim([6 11])

subplot(2,2,3)
plot(log10(eP(:)+1),log10(cP(:)+1),'.b')
title('P')
ylabel('FEISTY catch')
xlabel('Effort')

subplot(2,2,4)
plot(log10(eD(:)+1),log10(cD(:)+1),'.','color',[0 0.75 0.25])
title('D')
ylim([6 11])
print('-dpng',[ppath exper 'effort_scatter.png'])

%% Stats
%remove non-ocean grid cells
fid = find(~isnan(cPel));
eid = find(~isnan(eP));
pid = intersect(fid,eid);

fid = find(~isnan(cD));
eid = find(~isnan(eD));
did = intersect(fid,eid);

%r
[rP,pP]=corrcoef(log10(eP(pid)+1),log10(cPel(pid)+1));
[rD,pD]=corrcoef(log10(eD(did)+1),log10(cD(did)+1));
r(1) = rP(2,1);
r(2) = rD(2,1);

%root mean square error - Doesn't make sense with effort, not same units
o=log10(eP(pid)+1);
p=log10(cPel(pid)+1);
n = length(o);
num=nansum((p-o).^2);
rmse(1) = sqrt(num/n);

o=log10(eD(did)+1);
p=log10(cD(did)+1);
n = length(o);
num=nansum((p-o).^2);
rmse(2) = sqrt(num/n);

