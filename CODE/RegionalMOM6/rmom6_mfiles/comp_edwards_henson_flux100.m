% Compare Pippa Edwards POC flux at 100 m
% to Henson 2012 export flux

clear 
close all

%%
epath = '/Volumes/petrik-lab/Feisty/Obs_data/Export_stat_products/Pippa/';

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/OCB_BGC/export_flux/Edwards/';

%% load

EmgCd = EPCHn*1e3/365;

%% map info
[geolat,geolon] = meshgrid(lat,lon);

[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% maps

test = squeeze(meanlogpoc(:,:,1));

figure; pcolor(geolat); shading flat; colorbar
figure; pcolor(geolon); shading flat; colorbar
figure; pcolor(test); shading flat; colorbar

%%
close all
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,test)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-3 0]);
hcb = colorbar('h');
title('log10 mean POC flux at 100 m (mg/m2/day) Sep 1997')
stamp('')
print('-dpng',[ppath 'meanlogpoc_100m_test.png'])

%% mean over whole time
mpoc = mean(meanlogpoc,3,'omitnan');
yr=(double(date)/365)+1997+(9/12);

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,mpoc)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([2 5]);
hcb = colorbar('h');
title('log10 mean POC flux at 100 m (mg/m2/day) 1997-2022')
stamp('')
print('-dpng',[ppath 'meanlogpoc_100m_mean.png'])