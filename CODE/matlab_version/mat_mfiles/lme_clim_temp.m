% Calc LME mean temp
% Climatology
% Saved as mat files

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir = '/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([gpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])
load([cpath 'cobalt_zoop_biom_means.mat']);%,'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim')
load([cpath 'cobalt_det_biom_means.mat']);%,'det_mean_clim')

%%
ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

AREA_OCN = max(area,1);

%% Calc LMEs
tlme = lme_mask_onedeg;

lme_ptemp = NaN*ones(66,1);
lme_btemp = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %mean 
    lme_ptemp(L) = nanmean(ptemp_mean_clim(lid));
    lme_btemp(L) = nanmean(btemp_mean_clim(lid));
end

%%
save([gpath 'LME_clim_temp.mat'],'lme_ptemp','lme_btemp');

%% Figures
[ni,nj] = size(lon);
lme_pT = NaN*ones(ni,nj);
lme_bT = NaN*ones(ni,nj);

for L=1:66
    lid = find(tlme==L);

    lme_pT(lid) = lme_ptemp(L,1);
    lme_bT(lid) = lme_btemp(L,1);
    
end

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

%% 
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_pT)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 30]);
hcb = colorbar('h');
ylim(hcb,[-2 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean temp 0-100 m (^oC)')
print('-dpng',[gpath 'Clim_lme_Ptemp.png'])

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_bT)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 30]);
hcb = colorbar('h');
ylim(hcb,[-2 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean bottom temp (^oC)')
print('-dpng',[gpath 'Clim_lme_Btemp.png'])

