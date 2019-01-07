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
load([cpath 'cobalt_zoop_biom_means.mat'],'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim')
load([cpath 'cobalt_det_biom_means.mat'],'det_mean_clim')

%%
ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

z_mean = mz_mean_clim + lz_mean_clim;
z_loss = mzloss_mean_clim+lzloss_mean_clim;

AREA_OCN = max(area,1);

% mg/m2 --> total g
Az_mean = 1e-3 * z_mean .* AREA_OCN;
% mg/m2/d --> total g
Az_loss = 1e-3 * 365 * z_loss .* AREA_OCN;
% mg/m2/d --> total g; Mult flux by 365 - YES
Adet_mean = 1e-3 * 365 * det_mean_clim .* AREA_OCN;

%% Calc LMEs
tlme = lme_mask_onedeg;

lme_ptemp = NaN*ones(66,1);
lme_btemp = NaN*ones(66,1);
lme_z = NaN*ones(66,1);
lme_zl = NaN*ones(66,1);
lme_det = NaN*ones(66,1);
lme_az = NaN*ones(66,1);
lme_azl = NaN*ones(66,1);
lme_adet = NaN*ones(66,1);
lme_asz = NaN*ones(66,1);
lme_aszl = NaN*ones(66,1);
lme_asdet = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %mean 
    lme_ptemp(L) = nanmean(ptemp_mean_clim(lid));
    lme_btemp(L) = nanmean(btemp_mean_clim(lid));
    lme_z(L) = nanmean(z_mean(lid));
    lme_zl(L) = nanmean(z_loss(lid));
    lme_det(L) = nanmean(det_mean_clim(lid));
    lme_az(L) = nanmean(Az_mean(lid));
    lme_azl(L) = nanmean(Az_loss(lid));
    lme_adet(L) = nanmean(Adet_mean(lid));
    lme_asz(L) = nansum(Az_mean(lid));
    lme_aszl(L) = nansum(Az_loss(lid));
    lme_asdet(L) = nansum(Adet_mean(lid));
end

%%
save([gpath 'LME_clim_temp_zoop_det.mat'],'lme_ptemp','lme_btemp',...
    'lme_z','lme_zl','lme_det','lme_az','lme_azl','lme_adet',...
    'lme_asz','lme_aszl','lme_asdet');

%% Figures
[ni,nj] = size(lon);
lme_pT = NaN*ones(ni,nj);
lme_bT = NaN*ones(ni,nj);
lme_zp = NaN*ones(ni,nj);
lme_l = NaN*ones(ni,nj);
lme_d = NaN*ones(ni,nj);

for L=1:66
    lid = find(tlme==L);

    lme_pT(lid) = lme_ptemp(L,1);
    lme_bT(lid) = lme_btemp(L,1);
    lme_zp(lid) = lme_z(L,1);
    lme_l(lid) = lme_zl(L,1);
    lme_d(lid) = lme_det(L,1);
    
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
cmocean('thermal')
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
cmocean('thermal')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 30]);
hcb = colorbar('h');
ylim(hcb,[-2 30])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean bottom temp (^oC)')
print('-dpng',[gpath 'Clim_lme_Btemp.png'])
%%
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(lme_zp))
cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([2.5 3.5]);
hcb = colorbar('h');
ylim(hcb,[2.5 3.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean zoop biomass 0-100 m (log10 mgC m^-^2)')
print('-dpng',[gpath 'Clim_lme_zoop.png'])
%%
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(lme_l))
cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.5 2]);
hcb = colorbar('h');
ylim(hcb,[0.5 2])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean zoop loss 0-100 m (log10 mgC m^-^2 d^-^1)')
print('-dpng',[gpath 'Clim_lme_zloss.png'])

figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(lme_d))
cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.5 2]);
hcb = colorbar('h');
ylim(hcb,[0.5 2])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean bottom detritus flux (log10 mgC m^-^2 d^-^1)')
print('-dpng',[gpath 'Clim_lme_det.png'])

