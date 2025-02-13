% Get total NPP to calc TEs

clear all
close all

cpath='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
mos = [31,28,31,30,31,30,31,31,30,31,30,31]';

%% Climatology
%Fdet
load([cpath 'npp_100_1deg_ESM26_5yr_clim_191_195.mat'])

npp_100(npp_100<0) = 0;
npp_mean_clim=squeeze(nanmean(npp_100,1));

npp_tot = npp_100 .* repmat(mos,1,180,360);
npp_tot_clim=squeeze(nansum(npp_tot,1));

%%
save([cpath 'clim_npp_Dmeans_Ytot.mat'],'lon','lat',...
    'npp_mean_clim','npp_tot_clim');

%%
geolon_t = lon;
geolat_t = lat;

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_c,geolon_c,(npp_mean_clim))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 30]);
hcb = colorbar('h');
ylim(hcb,[-2 30])
set(gcf,'renderer','painters')
title('Top 100 m Temp Climatology (^oC)')
print('-dpng','t100_climatology.png')

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_c,geolon_c,(npp_tot_clim))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 30]);
hcb = colorbar('h');
ylim(hcb,[-2 30])
set(gcf,'renderer','painters')
title('Bottom Temp Climatology (^oC)')
print('-dpng','btemp_climatology.png')
