% P:D ratio by LME 
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));

%%
AREA_OCN = max(area,1);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%% Calc LMEs
tlme = lme_mask_onedeg;

load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_mbio','lme_area');

lme_area_km2 = lme_area * 1e-6;

% POEM LME biomass in MT
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

%% Figures

clme_mf = NaN*ones(180,360);
clme_mp = clme_mf;
clme_md = clme_mf;
clme_lp = clme_mf;
clme_ld = clme_mf;

plme_AllP = clme_mf;
plme_AllD = clme_mf;

lme_sf = NaN*ones(180,360);
lme_sp = lme_sf;
lme_sd = lme_sf;
lme_mf = lme_sf;
lme_mp = lme_sf;
lme_md = lme_sf;
lme_lp = lme_sf;
lme_ld = lme_sf;
lme_b = lme_sf;

for L=1:66
    lid = find(tlme==L);

    clme_mf(lid) = lme_mcatch(L,1);
    clme_mp(lid) = lme_mcatch(L,2);
    clme_md(lid) = lme_mcatch(L,3);
    clme_lp(lid) = lme_mcatch(L,4);
    clme_ld(lid) = lme_mcatch(L,5);

    lme_sf(lid) = lme_mbio(L,1);
    lme_sp(lid) = lme_mbio(L,2);
    lme_sd(lid) = lme_mbio(L,3);
    lme_mf(lid) = lme_mbio(L,4);
    lme_mp(lid) = lme_mbio(L,5);
    lme_md(lid) = lme_mbio(L,6);
    lme_lp(lid) = lme_mbio(L,7);
    lme_ld(lid) = lme_mbio(L,8);
    lme_b(lid) = lme_mbio(L,9);
    
    plme_AllP(lid) = plme_Pmcatch(L);
    plme_AllD(lid) = plme_Dmcatch(L);
end

clme_AllP = clme_mp+clme_lp;
clme_AllD = clme_md+clme_ld;

lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;
lme_AllM = lme_mf+lme_mp+lme_md;
lme_AllL = lme_lp+lme_ld;

%% Ratios
rPD_biom = lme_AllP ./ (lme_AllP+lme_AllD);
rPF_biom = lme_AllP ./ (lme_AllP+lme_AllF);
rLM_biom = lme_AllL ./ (lme_AllL+lme_AllM);
rPD_catch = clme_AllP ./ (clme_AllP+clme_AllD);
rPD_catch_mtkm2 = plme_AllP ./ (plme_AllP+plme_AllD);

save([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'rPD_biom','rPF_biom','rLM_biom','rPD_catch','rPD_catch_mtkm2',...
    '-append');

%% Plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

revamp=colormap(cmocean('amp'));
revamp=flipud(revamp);

%% Biomass

figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPD_biom)
cmocean('balance')
%cmocean('amp')
%colormap(revamp)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
hcb = colorbar('h');
ylim(hcb,[0 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean ratio P:D biomass (g)')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_ratioPD_biom.png'])

%% Catch
% grams
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPD_catch)
cmocean('balance')
%cmocean('amp')
%colormap(revamp)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
hcb = colorbar('h');
ylim(hcb,[0 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean ratio P:D catch (g)')
stamp(cfile)
%print('-dpng',[ppath 'Clim_fished_',harv,'_LME_ratioPD_catch.png'])

% MT/km2
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPD_catch)
cmocean('balance')
%cmocean('amp')
%colormap(revamp)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
hcb = colorbar('h');
ylim(hcb,[0 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology LME mean ratio P:D catch (MT km^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_ratioPD_catch_mtkm2.png'])

%% 3 figure subplot P:D, P:F, M:L
figure(30)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPD_biom)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rPF_biom)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,rLM_biom)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Clim_fished_' harv '_LME_ratios_subplot.png'])

