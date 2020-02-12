% Visualize difference between
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100
% Biome size and biomass in each biome

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([bpath 'COBALT_biomes_50yr_means_5100.mat'],'biome_hist','biome_fore');

load([cpath 'hindcast_gridspec.mat'],'AREA_OCN','geolon_t','geolat_t');
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);

% colors
cmPG=cbrewer('div','PiYG',50,'PCHIP');
cmPO=cbrewer('div','PuOr',50,'PCHIP');
cmG=cbrewer('seq','Greens',50,'PCHIP');
cmB=cbrewer('seq','Blues',50,'PCHIP');
cmP=cbrewer('seq','Purples',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmBW=cbrewer('seq','Greys',50,'PCHIP');

cmapD = cmR(35,:);
cmapD(2,:) = cmB(35,:);
cmapD(3,:) = cmG(35,:);

%% Hindcast grid
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

% Hindcast
load([fpath 'Biomes_Hist_',harv,'_' cfile '.mat']);

% Forecast
load([fpath 'Biomes_Fore_',harv,'_' cfile '.mat']);

Hbiome_mbF = Hbiome_mbio50(:,1) + Hbiome_mbio50(:,4);
Hbiome_mbP = Hbiome_mbio50(:,2) + Hbiome_mbio50(:,5) + Hbiome_mbio50(:,7);
Hbiome_mbD = Hbiome_mbio50(:,3) + Hbiome_mbio50(:,6) + Hbiome_mbio50(:,8);
Hbiome_mbB = Hbiome_mbio50(:,9);
Hbiome_mbA = Hbiome_mbF + Hbiome_mbP + Hbiome_mbD;

Fbiome_mbF = Fbiome_mbio50(:,1) + Fbiome_mbio50(:,4);
Fbiome_mbP = Fbiome_mbio50(:,2) + Fbiome_mbio50(:,5) + Fbiome_mbio50(:,7);
Fbiome_mbD = Fbiome_mbio50(:,3) + Fbiome_mbio50(:,6) + Fbiome_mbio50(:,8);
Fbiome_mbB = Fbiome_mbio50(:,9);
Fbiome_mbA = Fbiome_mbF + Fbiome_mbP + Fbiome_mbD;

FbiomeH_mbF = FbiomeH_mbio50(:,1) + FbiomeH_mbio50(:,4);
FbiomeH_mbP = FbiomeH_mbio50(:,2) + FbiomeH_mbio50(:,5) + FbiomeH_mbio50(:,7);
FbiomeH_mbD = FbiomeH_mbio50(:,3) + FbiomeH_mbio50(:,6) + FbiomeH_mbio50(:,8);
FbiomeH_mbB = FbiomeH_mbio50(:,9);
FbiomeH_mbA = FbiomeH_mbF + FbiomeH_mbP + FbiomeH_mbD;

%% Map grid

Hbiome_AllF = NaN*ones(360,200);
Hbiome_AllP = Hbiome_AllF;
Hbiome_AllD = Hbiome_AllF;
Hbiome_AllB = Hbiome_AllF;
Hbiome_All  = Hbiome_AllF;

FbiomeH_AllF = NaN*ones(360,200);
FbiomeH_AllP = FbiomeH_AllF;
FbiomeH_AllD = FbiomeH_AllF;
FbiomeH_AllB = FbiomeH_AllF;
FbiomeH_All  = FbiomeH_AllF;

Fbiome_AllF = NaN*ones(360,200);
Fbiome_AllP = Fbiome_AllF;
Fbiome_AllD = Fbiome_AllF;
Fbiome_AllB = Fbiome_AllF;
Fbiome_All  = Fbiome_AllF;

for L=1:3
    hid = find(biome_hist==L);
    fid = find(biome_fore==L);
    
    Hbiome_AllF(hid) = Hbiome_mbF(L);
    Hbiome_AllP(hid) = Hbiome_mbP(L);
    Hbiome_AllD(hid) = Hbiome_mbD(L);
    Hbiome_AllB(hid) = Hbiome_mbB(L);
    Hbiome_All(hid)  = Hbiome_mbA(L);
    
    FbiomeH_AllF(hid) = FbiomeH_mbF(L);
    FbiomeH_AllP(hid) = FbiomeH_mbP(L);
    FbiomeH_AllD(hid) = FbiomeH_mbD(L);
    FbiomeH_AllB(hid) = FbiomeH_mbB(L);
    FbiomeH_All(hid)  = FbiomeH_mbA(L);
    
    Fbiome_AllF(fid) = Fbiome_mbF(L);
    Fbiome_AllP(fid) = Fbiome_mbP(L);
    Fbiome_AllD(fid) = Fbiome_mbD(L);
    Fbiome_AllB(fid) = Fbiome_mbB(L);
    Fbiome_All(fid)  = Fbiome_mbA(L);
end


%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

%%
pdiffF = (Fbiome_AllF-Hbiome_AllF) ./ Hbiome_AllF;
pdiffP = (Fbiome_AllP-Hbiome_AllP) ./ Hbiome_AllP;
pdiffD = (Fbiome_AllD-Hbiome_AllD) ./ Hbiome_AllD;
pdiffAll = (Fbiome_All-Hbiome_All) ./ Hbiome_All;

diffF = (Fbiome_AllF-Hbiome_AllF) ;
diffP = (Fbiome_AllP-Hbiome_AllP) ;
diffD = (Fbiome_AllD-Hbiome_AllD) ;
diffAll = (Fbiome_All-Hbiome_All) ;

HpdiffF = (FbiomeH_AllF-Hbiome_AllF) ./ Hbiome_AllF;
HpdiffP = (FbiomeH_AllP-Hbiome_AllP) ./ Hbiome_AllP;
HpdiffD = (FbiomeH_AllD-Hbiome_AllD) ./ Hbiome_AllD;
HpdiffAll = (FbiomeH_All-Hbiome_All) ./ Hbiome_All;

HdiffF = (FbiomeH_AllF-Hbiome_AllF) ;
HdiffP = (FbiomeH_AllP-Hbiome_AllP) ;
HdiffD = (FbiomeH_AllD-Hbiome_AllD) ;
HdiffAll = (FbiomeH_All-Hbiome_All) ;

%% F
figure(1)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Hbiome_AllF))
colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Hindcast F');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Fbiome_AllF))
colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Forecast F');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('% Forecast - Hindcast F');
print('-dpng',[pp 'Hist_Fore_' harv '_biome_pdiffF_cmR.png'])

%% P
figure(2)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Hbiome_AllP))
colormap(cmB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Hindcast P');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Fbiome_AllP))
colormap(cmB)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Forecast P');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('% Forecast - Hindcast P');
print('-dpng',[pp 'Hist_Fore_' harv '_biome_pdiffP_cmB.png'])

%% D
figure(3)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Hbiome_AllD))
colormap(cmG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Hindcast D');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Fbiome_AllD))
colormap(cmG)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Forecast D');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('% Forecast - Hindcast D');
print('-dpng',[pp 'Hist_Fore_' harv '_biome_pdiffP_cmG.png'])

%% All
figure(4)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Hbiome_All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Hindcast All');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Fbiome_All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('log_1_0 Forecast All');

% Diff
subplot('Position',[0.5 0.3 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('% Forecast - Hindcast All');
print('-dpng',[pp 'Hist_Fore_' harv '_biome_pdiffP_cmB.png'])

%% All 4 on subplots 
figure(5)
% F
subplot('Position',[0.01 0.52 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast F');

% P
subplot('Position',[0.51 0.52 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast P');
colorbar('Position',[0.25 0.5 0.5 0.03],'orientation','horizontal')

% D
subplot('Position',[0.01 0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast D');

% All
subplot('Position',[0.51 0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,pdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast All');

print('-dpng',[pp 'Hist_Fore_' harv '_biome_pdiff_subplot.png'])

%% All 4 on subplots using historic biome areas
figure(6)
% F
subplot('Position',[0.01 0.52 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,HpdiffF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast F');

% P
subplot('Position',[0.51 0.52 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,HpdiffP)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast P');
colorbar('Position',[0.25 0.5 0.5 0.03],'orientation','horizontal')

% D
subplot('Position',[0.01 0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,HpdiffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast D');

% All
subplot('Position',[0.51 0 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,HpdiffAll)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
set(gcf,'renderer','painters')
title('% Forecast - Hindcast All');

print('-dpng',[pp 'Hist_Fore_' harv '_biome_hist_pdiff_subplot.png'])


%% regions with change in biome def
figure(7)
subplot('Position',[0 0.53 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,biome_hist)
colormap(cmapD)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 3]);
%colorbar('southoutside','Ticks',[1 2 3],'TickLabels',{'F','P','D'})
colorbar('Position',[0.05 0.53 0.4 0.025],'orientation','horizontal',...
    'Ticks',[1 2 3],'TickLabels',{'LC','ECCS','ECSS'})
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Hindcast');

subplot('Position',[0 0.09 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,biome_fore)
colormap(cmapD)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 3]);
% colorbar
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
title('Forecast');

% Diff
subplot('Position',[0.5 0.33 0.475 0.33])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(geolat_t,geolon_t,biome_fore-biome_hist)
cmocean('balance')
%colormap(cmR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
colorbar('Position',[0.95 0.34 0.02 0.33],'orientation','vertical','Ticks',[-1 0 1])
set(gcf,'renderer','painters')
set(gca,'XColor','none','YColor','none')
%text(-0.1,1.5,'Large pelagics','FontWeight','bold','HorizontalAlignment','center');

print('-dpng',[pp 'Hist_Fore_',harv,'_biome_change.png'])

