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
cfile = 'Dc_enc70-b200_m4-b175-k08_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_nu_gam_die_clev_Climatol_' harv '_' cfile '.mat']);

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

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

Zsf(ID)=sf_mdie;
Zsp(ID)=sp_mdie;
Zsd(ID)=sd_mdie;
Zmf(ID)=mf_mdie;
Zmp(ID)=mp_mdie;
Zmd(ID)=md_mdie;
Zlp(ID)=lp_mdie;
Zld(ID)=ld_mdie;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mdie));

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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_SP_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_SF_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_SD_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_MP_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_MF_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_MD_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_LP_die.png'])
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
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_LD_die.png'])

%% Diff maps of all fish
All = (Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld)/8;
AllF = (Zsf+Zmf)/2;
AllP = (Zsp+Zmp+Zlp)/3;
AllD = (Zsd+Zmd+Zld)/3;
AllS = (Zsp+Zsf+Zsd)/3;
AllM = (Zmp+Zmf+Zmd)/3;
AllL = (Zlp+Zld)/2;

%% ALL
% figure(21)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(All))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All fishes')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_All_die.png'])
% 
% % all F
% figure(22)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(AllF))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All F')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllF_die.png'])
% 
% % all D
% figure(23)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(AllD))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All D')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllD_die.png'])
% 
% % All P
% figure(24)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(AllP))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All P')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllP_die.png'])

% All 4 on subplots
figure(27)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title(' mean All F')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
set(gcf,'renderer','painters')
title(' mean All D')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
set(gcf,'renderer','painters')
title(' mean All P')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
set(gcf,'renderer','painters')
title(' mean All fishes')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot_die.png'])

%%
% % all S
% figure(2)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(AllS))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All S')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllS_die.png'])
% 
% % all M
% figure(3)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(AllM))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All M')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllM_die.png'])
% 
% % All L
% figure(4)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,(AllL))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0.2 1]);
% hcb = colorbar('h');
% ylim(hcb,[0.2 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology  mean All L')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_global_AllL_die.png'])

% All 4 on subplots
figure(5)
% all S
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(AllS))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title(' mean All S')

% all L
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(AllL))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
set(gcf,'renderer','painters')
title(' mean All L')

% All M
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(AllM))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
set(gcf,'renderer','painters')
title(' mean All M')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 0.025]);
set(gcf,'renderer','painters')
title(' mean All fishes')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_sizes_subplot_die.png'])


