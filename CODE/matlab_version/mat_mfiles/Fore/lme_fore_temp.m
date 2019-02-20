% Calc LME mean temp
% Forecast 2006-2100
% Last 50 years?
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
load([cpath 'cobalt_temp_means.mat']);

%% 
ptemp_mean_hist=ptemp_mean_fore;
btemp_mean_hist=btemp_mean_fore;

ID = grid(:,1);

%% Calc LMEs
tlme = lme_mask_esm2m';

lme_ptemp = NaN*ones(66,1);
lme_btemp = NaN*ones(66,1);


for L=1:66
    lid = find(tlme==L);
    %mean 
    lme_ptemp(L) = nanmean(ptemp_mean_hist(lid));
    lme_btemp(L) = nanmean(btemp_mean_hist(lid));
    
end

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
    
end

%%
save([cpath 'LME_fore_temp.mat'],'lme_ptemp','lme_btemp','lme_pT','lme_bT');


%% plot info
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
title('Forecast 2051-2100 LME mean temp 0-100 m (^oC)')
print('-dpng',[gpath 'Fore_lme_Ptemp.png'])

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
title('Forecast 2051-2100 LME mean bottom temp (^oC)')
print('-dpng',[gpath 'Fore_lme_Btemp.png'])

