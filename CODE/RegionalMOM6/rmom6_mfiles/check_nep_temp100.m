% Check TP calc with calc thkcello

clear
close all

cpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';

%% Units
%poc flux: mol N m-2 s-1
%zoo: mol N m-2
%tp: degC
%tb: degC

load([cpath 'temp_100.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat'],'temp_100');

%%
temp_100(temp_100 > 1.0e19) = nan;

%%
load([cpath 'nep_raw_ocean_static_gridspec.mat'],'geolon','geolat');

[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);

%% NEPac
plotminlat=10; 
plotmaxlat=85;
plotminlon=156;
plotmaxlon=-104;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%%
figure
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,squeeze(temp_100(:,:,1)))
cmocean('thermal')
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30]);
hcb = colorbar('h');
title('t100')