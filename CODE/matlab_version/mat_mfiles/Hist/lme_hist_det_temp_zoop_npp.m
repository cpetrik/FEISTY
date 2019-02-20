% Calc LME mean temp
% Historic 1951-2000
% Saved as mat files

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm2m_hist/';
cdir = '/Volumes/GFDL/GCM_DATA/ESM2M_hist/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
load([gpath 'lme_mask_esm2m.mat']);
grid = csvread([gpath 'grid_csv.csv']);
load([cpath 'cobalt_det_biom_means.mat']);
load([cpath 'cobalt_npp_means.mat']);
load([cpath 'cobalt_temp_means.mat']);
load([cpath 'cobalt_zoop_biom_means.mat']);

%% CHECK UNITS
ptemp_mean_hist=ptemp_mean_hist;
btemp_mean_hist=btemp_mean_hist;

% molN/m2 --> g/m2
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mz_mean_hist = mz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
lz_mean_hist = lz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
% molN/m2/s --> g/m2/d
mzloss_mean_hist = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_mean_hist = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_mean_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_mean_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

z_mean = mz_mean_hist + lz_mean_hist;
z_loss = mzloss_mean_hist + lzloss_mean_hist;

%AREA_OCN = max(AREA_OCN,1); Not sure what units area is in ~10^-5 vs. 10^9
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

ID = grid(:,1);

% g/m2 --> total g
Az_mean = z_mean .* AREA_OCN;
% g/m2/d --> total g; Mult flux by 365 - YES
Az_loss = 365 * z_loss .* AREA_OCN;
Adet_mean = 365 * det_mean_hist .* AREA_OCN;
Anpp_mean = 365 * npp_mean_hist .* AREA_OCN;

%% Calc LMEs
tlme = lme_mask_esm2m';

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
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %mean 
    lme_ptemp(L) = nanmean(ptemp_mean_hist(lid));
    lme_btemp(L) = nanmean(btemp_mean_hist(lid));
    lme_z(L) = nanmean(z_mean(lid));
    lme_zl(L) = nanmean(z_loss(lid));
    lme_det(L) = nanmean(det_mean_hist(lid));
    lme_az(L) = nanmean(Az_mean(lid));
    lme_azl(L) = nanmean(Az_loss(lid));
    lme_adet(L) = nanmean(Adet_mean(lid));
    lme_asz(L) = nansum(Az_mean(lid));
    lme_aszl(L) = nansum(Az_loss(lid));
    lme_asdet(L) = nansum(Adet_mean(lid));
    lme_area(L) = nansum(AREA_OCN(lid));
end

%%
save([cpath 'LME_hist_temp_zoop_det.mat'],'lme_ptemp','lme_btemp',...
    'lme_z','lme_zl','lme_det','lme_az','lme_azl','lme_adet',...
    'lme_asz','lme_aszl','lme_asdet','lme_area');

%% Figures
[ni,nj] = size(geolon_t);
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
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
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
title('Historic 1951-2000 LME mean temp 0-100 m (^oC)')
print('-dpng',[cpath 'Hist_lme_Ptemp.png'])

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
title('Historic 1951-2000 LME mean bottom temp (^oC)')
print('-dpng',[cpath 'Hist_lme_Btemp.png'])

%%
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(lme_zp))
cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.25 1.25]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Historic 1951-2000 LME mean zoop biomass 0-100 m (log10 g m^-^2)')
print('-dpng',[cpath 'Hist_lme_zoop.png'])

%%
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(lme_l))
cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Historic 1951-2000 LME mean zoop loss 0-100 m (log10 g m^-^2 d^-^1)')
print('-dpng',[cpath 'Hist_lme_zloss.png'])

figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(lme_d))
cmocean('matter')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Historic 1951-2000 LME mean bottom detritus flux (log10 g m^-^2 d^-^1)')
print('-dpng',[cpath 'Hist_lme_det.png'])

