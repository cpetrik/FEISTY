% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
%harv = 'fish_F030_P060_D060';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

%load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'Climatol_All_fish03_ICbiom_1e-10.mat']);
%icbt = {'1e-10';'1e-7';'1e-6';'1e-4';'1e-3';'1e-2';'1e3'};
icbt = '1e-10';

close all

% plot info
[ni,nj]=size(lon);
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

geolat_t=lat;
geolon_t=lon;

% colors
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;


ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mean));

% bent
figure(50)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.5 0.5]);
hcb = colorbar('h');
ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 mean Benthic inverts (g m^-^2)')
stamp(['ICbiom_' icbt '_' harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_BENT_ICbiom_',icbt,'.png'])

%
% mgZb = (Zb/9)*1e3;
% figure(51)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(mgZb))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-0.8 2.3]);
% hcb = colorbar('h');
% ylim(hcb,[-0.8 2.3])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean Benthic inverts (mg m^-^2)')
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_BENT_mgC_ICbiom_',icbt,'.png'])

%% sp
% figure(11)
% surf(lon,lat,log10(Zsp)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Larval P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_SP_ICbiom_',icbt,'.png'])
%
% % sf
% figure(12)
% surf(lon,lat,log10(Zsf)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Larval F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_SF_ICbiom_',icbt,'.png'])
%
% % sd
% figure(13)
% surf(lon,lat,log10(Zsd)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Larval D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_SD_ICbiom_',icbt,'.png'])
%
% % mp
% figure(14)
% surf(lon,lat,log10(Zmp)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Juvenile P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_MP_ICbiom_',icbt,'.png'])
%
% % mf
% figure(15)
% surf(lon,lat,log10(Zmf)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Adult F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_MF_ICbiom_',icbt,'.png'])
%
% % md
% figure(16)
% surf(lon,lat,log10(Zmd)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Juvenile D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_MD_ICbiom_',icbt,'.png'])
%
% % lp
% figure(17)
% surf(lon,lat,log10(Zlp)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Adult P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_LP_ICbiom_',icbt,'.png'])
%
% % ld
% figure(18)
% surf(lon,lat,log10(Zld)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Adult D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_LD_ICbiom_',icbt,'.png'])

% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllM+AllL);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
FracPDs = Zsp ./ (Zsp+Zsd);
FracPDm = Zmp ./ (Zmp+Zmd);
FracPDl = Zlp ./ (Zlp+Zld);
FracPFs = Zsp ./ (Zsp+Zsf);
FracPFm = Zmp ./ (Zmp+Zmf);
FracPFvDs = (Zsp+Zsf) ./ (Zsp+Zsf+Zsd);
FracPFvDm = (Zmp+Zmf) ./ (Zmp+Zmf+Zmd);

% ALL
% figure(21)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(All))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 2]);
% hcb = colorbar('h');
% ylim(hcb,[-1 2])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All fishes (g m^-^2)')
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_All_ICbiom_',icbt,'.png'])
% 
% % all F
% figure(22)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllF))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% hcb = colorbar('h');
% ylim(hcb,[-1 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All F (g m^-^2)')
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllF_ICbiom_',icbt,'.png'])
% 
% % all D
% figure(23)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllD))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% hcb = colorbar('h');
% ylim(hcb,[-1 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All D (g m^-^2)')
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllD_ICbiom_',icbt,'.png'])
% 
% % All P
% figure(24)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllP))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% hcb = colorbar('h');
% ylim(hcb,[-1 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All P (g m^-^2)')
% stamp(['ICbiom_' icbt '_' harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllP_ICbiom_',icbt,'.png'])

%% All 4 on subplots
figure(27)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[-1 1])                   %Set color axis if needed
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[-1 1])
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
%     stamp(['ICbiom_' icbt '_' harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot_ICbiom_',icbt,'.png'])


%% 3 figure subplot P:D, P:F, M:L
figure(30)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
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
surfm(geolat_t,geolon_t,FracPF)
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
surfm(geolat_t,geolon_t,FracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
%stamp(['ICbiom_' icbt '_' harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_ratios_subplot_v3_ICbiom_',icbt,'.png'])


