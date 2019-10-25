% P:D ratio by US LMEs - no Arctic
% Climatology
% 150 years
% Saved as mat files
% Compare to SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);

tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

%% FEISTY 
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/US/'];
dpath = [dp cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_area','rPD_biom','rPD_catch','rPD_catch_mtkm2');

lme_area_km2 = lme_area * 1e-6;

%% FEISTY  LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch);
l10pF=log10(plme_Fmcatch);
l10pP=log10(plme_Pmcatch);
l10pD=log10(plme_Dmcatch);

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = [1:7,10,65];

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% Put both on grid
tlme = lme_mask_onedeg;

pF_grid = NaN*ones(180,360);
pP_grid = NaN*ones(180,360);
pD_grid = NaN*ones(180,360);
pAll_grid = NaN*ones(180,360);
sF_grid = NaN*ones(180,360);
sP_grid = NaN*ones(180,360);
sD_grid = NaN*ones(180,360);
sAll_grid = NaN*ones(180,360);
pFracPD_grid = NaN*ones(180,360);
sFracPD_grid = NaN*ones(180,360);

for n=1:length(keep)
    L=keep(n);
    lid = find(tlme==L);
    pF_grid(lid) = l10pF(L);
    pP_grid(lid) = l10pP(L);
    pD_grid(lid) = l10pD(L);
    pAll_grid(lid) = l10p(L);
    pFracPD_grid(lid) = pFracPD(L);
    
    sF_grid(lid) = l10sF(L);
    sP_grid(lid) = l10sP(L);
    sD_grid(lid) = l10sD(L);
    sAll_grid(lid) = l10s(L);
    sFracPD_grid(lid) = sFracPD(L);
end

%% Comparison 

diffF = pF_grid - sF_grid;
diffP = pP_grid - sP_grid;
diffD = pD_grid - sD_grid;
diffAll = pAll_grid - sAll_grid;
diffPD = pFracPD_grid - sFracPD_grid;

%% Plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=10; %Set these bounds for your data
plotmaxlat=65;
plotminlon=-180;
plotmaxlon=-50;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

cmYOR=cbrewer('seq','YlOrRd',28);
cmRP=cbrewer('seq','RdPu',28);
cmPR=cbrewer('seq','PuRd',28);

x=-6:0.1:3;

%% Subplot with maps and corr no temp color
% Forage
figure(1)
% FEISTY
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,pF_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3.5 1.5]);
colorbar('Position',[0.25 0.575 0.5 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'A','FontSize',14)

% SAUP
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,sF_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3.5 1.5]);
%colorbar('Position',[0.94 0.575 0.5 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'B','FontSize',14)

% Diff
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.075 0.45 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'C','FontSize',14)

%Corr
subplot('Position',[0.58 0.1 0.38 0.38])
plot(x,x,'--k');hold on;
scatter(l10sF(keep),l10pF(keep),50,'k','filled'); hold on;
axis([-6 1 -6 1])
xlabel('SAU log_1_0 F catch (MT)')
ylabel('FEISTY log_1_0 F catch (MT)')
text(-5.8,0.5,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_US_LMEs_noArc_F_SAUP_comp_subplot_corr.png'])

%% Large Pel
figure(2)
% FEISTY
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,pP_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
colorbar('Position',[0.25 0.575 0.5 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'A','FontSize',14)

% SAUP
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,sP_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'B','FontSize',14)

% Diff
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.075 0.45 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'C','FontSize',14)

%Corr
subplot('Position',[0.58 0.1 0.38 0.38])
plot(x,x,'--k');hold on;
scatter(l10sP(keep),l10pP(keep),50,'k','filled'); hold on;
axis([-3 1 -3 1])
xlabel('SAU log_1_0 P catch (MT)')
ylabel('FEISTY log_1_0 P catch (MT)')
text(-2.75,0.5,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_US_LMEs_noArc_P_SAUP_comp_subplot_corr.png'])

%% Demersal
figure(3)
% FEISTY
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,pD_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
colorbar('Position',[0.25 0.575 0.5 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'A','FontSize',14)

% SAUP
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,sD_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'B','FontSize',14)

% Diff
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.075 0.45 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'C','FontSize',14)

%Corr
subplot('Position',[0.58 0.1 0.38 0.38])
plot(x,x,'--k');hold on;
scatter(l10sD(keep),l10pD(keep),50,'k','filled'); hold on;
axis([-1.5 1.5 -1.5 1.5])
xlabel('SAU log_1_0 D catch (MT)')
ylabel('FEISTY log_1_0 D catch (MT)')
text(-1.3,1.2,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_US_LMEs_noArc_D_SAUP_comp_subplot_corr.png'])

%% All
figure(4)
% FEISTY
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,pAll_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1.5]);
colorbar('Position',[0.25 0.575 0.5 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'A','FontSize',14)

% SAUP
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,sAll_grid)
colormap(cmYOR);
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1.5]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'B','FontSize',14)

% Diff
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.075 0.45 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'C','FontSize',14)

%Corr
subplot('Position',[0.58 0.1 0.38 0.38])
plot(x,x,'--k');hold on;
scatter(l10s(keep),l10p(keep),50,'k','filled'); hold on;
axis([-2 1.5 -2 1.5])
xlabel('SAU log_1_0 All catch (MT)')
ylabel('FEISTY log_1_0 All catch (MT)')
text(-1.75,1.2,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_US_LMEs_noArc_All_SAUP_comp_subplot_corr.png'])

%% P:D
figure(5)
% FEISTY
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,pFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.575 0.5 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'A','FontSize',14)

% SAUP
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,sFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'B','FontSize',14)

% Diff
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.025 0.075 0.45 0.025],'orientation','horizontal')                   
set(gcf,'renderer','painters')
text(-1.3,1.4,'C','FontSize',14)

%Corr
subplot('Position',[0.58 0.1 0.38 0.38])
plot(x,x,'--k');hold on;
scatter(sFracPD(keep),pFracPD(keep),50,'k','filled'); hold on;
axis([0 1 0 1])
xlabel('SAU fraction large pelagics in catch')
ylabel('FEISTY fraction large pelagics in catch')
text(0.05,0.9,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_US_LMEs_noArc_PD_SAUP_comp_subplot_corr.png'])



%% Subplot with all diff maps 
figure(6)
%F
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'A','FontSize',14)

% P
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'B','FontSize',14)

% D
subplot('Position',[0.0 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'C','FontSize',14)

%All
subplot('Position',[0.5 0.03 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,diffAll)
cmocean('balance')
colorbar('Position',[0.25 0.525 0.5 0.025],'orientation','horizontal')                   
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
text(-1.3,1.4,'D','FontSize',14)
print('-dpng',[ppath 'Clim_' harv '_US_LMEs_noArc_AllDiffs_SAUP_comp_subplot.png'])

