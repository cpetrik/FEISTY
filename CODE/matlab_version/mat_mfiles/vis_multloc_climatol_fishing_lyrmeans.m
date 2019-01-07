% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 100 years
% Only last 12 months saved to compare to SAUP for tuning

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
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

%% colors
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

%%
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE10_CC100_lgRE00100_mdRE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3';

fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Clim_means_' harv '.mat']);

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

close all

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

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;

mf_my(mf_my<0)=0;
Cmf(ID)=mf_my;
Cmp(ID)=mp_my;
Cmd(ID)=md_my;
Clp(ID)=lp_my;
Cld(ID)=ld_my;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mean));

%% ocean cells
% figure(55)
% surf(lon,lat,ocean); view(2); hold on;
% shading flat
% title('Water cells')
% colormap('jet')
% colorbar('h')
% caxis([1 2])
% stamp(cfile)
% print('-dpng',[ppath 'Ocean_cells.png'])

% bent
% figure(50)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(Zb))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2.5 0.5]);
% hcb = colorbar('h');
% ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean Benthic inverts (g m^-^2)')
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_BENT.png'])
% 
% %
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_BENT_mgC.png'])

% % sp
% figure(11)
% surf(lon,lat,log10(Zsp)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Larval P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_SP.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_SF.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_SD.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_MP.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_MF.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_MD.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_LP.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_LD.png'])

%% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
FracPDs = Zsp ./ (Zsp+Zsd);
FracPDm = Zmp ./ (Zmp+Zmd);
FracPDl = Zlp ./ (Zlp+Zld);
FracPFs = Zsp ./ (Zsp+Zsf);
FracPFm = Zmp ./ (Zmp+Zmf);
FracPFvDs = (Zsp+Zsf) ./ (Zsp+Zsf+Zsd);
FracPFvDm = (Zmp+Zmf) ./ (Zmp+Zmf+Zmd);

%% ALL
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_All.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_AllF.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_AllD.png'])
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
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_AllP.png'])
% 
% % FracPD
% figure(25)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,FracPD)
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 1]);
% hcb = colorbar('h');
% ylim(hcb,[0 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology Fraction P:D')
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPD.png'])
% 
% % FracPF
% figure(26)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,FracPF)
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([0 1]);
% hcb = colorbar('h');
% ylim(hcb,[0 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology Fraction P:F')
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPF.png'])

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
%     stamp(cfile)
print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot.png'])


%% cmocean
% % All 4 on subplots thermal
% figure
% % all F
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllF))
% cmocean('thermal')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('A. Forage Fishes')
% 
% % all D
% subplot('Position',[0 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllD))
% cmocean('thermal')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('C. Demersals')
% 
% % All P
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllP))
% cmocean('thermal')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('B. Large Pelagics')
% 
% % All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(All))
% cmocean('thermal')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('D. All fishes')
% print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot_thermal.png'])
% 
% %% All 4 on subplots matter
% figure
% % all F
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllF))
% cmocean('matter')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('A. Forage Fishes')
% 
% % all D
% subplot('Position',[0 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllD))
% cmocean('matter')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('C. Demersals')
% 
% % All P
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllP))
% cmocean('matter')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('B. Large Pelagics')
% 
% % All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(All))
% cmocean('matter')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('D. All fishes')
% print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot_matter.png'])
% 
% 
% %% All 4 on subplots haline
% figure
% % all F
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllF))
% cmocean('haline')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('A. Forage Fishes')
% 
% % all D
% subplot('Position',[0 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllD))
% cmocean('haline')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('C. Demersals')
% 
% % All P
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(AllP))
% cmocean('haline')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('B. Large Pelagics')
% 
% % All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,log10(All))
% cmocean('haline')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% set(gcf,'renderer','painters')
% title('D. All fishes')
% print('-dpng',[ppath 'Climatol_' harv '_global_All_subplot_haline.png'])

%% Ratios on subplots
figure(28)
% all P:F
subplot('Position',[0 0.55 1 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPF)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[0 1])
colorbar('Position',[0.2 0.475 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('A. Large Pelagics : Forage Fishes')

% all P:D
subplot('Position',[0 0 1 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('B. Large Pelagics : Demersals')
%     stamp(cfile)
print('-dpng',[ppath 'Climatol_' harv '_global_ratios_subplot.png'])

%% Ratios on subplots red-white-blue
figure(29)
% all P:F
subplot('Position',[0 0.55 1 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPF)
%     colormap(cmap_color_rb)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[0 1])
colorbar('Position',[0.2 0.475 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('A. Large Pelagics : Forage Fishes')

% all P:D
subplot('Position',[0 0 1 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
%     colormap(cmap_color_rb)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('B. Large Pelagics : Demersals')
%     stamp(cfile)
print('-dpng',[ppath 'Climatol_' harv '_global_ratios_subplot_v2.png'])


%% CATCH
% mp
% figure(34)
% surf(lon,lat,log10(Cmp)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Juvenile P catch (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-5 -2])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_MP_catch.png'])

% mf
%         figure(35)
%         surf(lon,lat,log10(Cmf)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean Adult F catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_MF_catch.png'])
%
% md
% figure(36)
% surf(lon,lat,log10(Cmd)); view(2); hold on;
% shading flat
% xlim([0 360])
% ylim([-90 90])
% title('log10 mean Juvenile D catch (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-5 -2])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_MD_catch.png'])

% lp
%         figure(37)
%         surf(lon,lat,log10(Clp)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean Adult P catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_LP_catch.png'])
%
%         % ld
%         figure(38)
%         surf(lon,lat,log10(Cld)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean Adult D catch (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_LD_catch.png'])
%
%         %% Diff maps of all fish
%         CAll = Cmp+Cmf+Cmd+Clp+Cld;
%         CAllF = Cmf;
%         CAllP = Cmp+Clp;
%         CAllD = Cmd+Cld;
%         CAllM = Cmp+Cmf+Cmd;
%         CAllL = Clp+Cld;
%
%         % ALL
%         figure(39)
%         surf(lon,lat,log10(CAll)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean catch All Fishes (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_All_catch.png'])
%
%         % all F
%         figure(40)
%         surf(lon,lat,log10(CAllF)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean catch All F (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_AllF_catch.png'])
%
%         % all D
%         figure(41)
%         surf(lon,lat,log10(CAllD)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean catch All D (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_AllD_catch.png'])
%
%         % All P
%         figure(42)
%         surf(lon,lat,log10(CAllP)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean catch All P (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_AllP_catch.png'])
%
%         % all M
%         figure(43)
%         surf(lon,lat,log10(CAllM)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean catch All M (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_AllM_catch.png'])
%
%         % All L
%         figure(44)
%         surf(lon,lat,log10(CAllL)); view(2); hold on;
%         shading flat
%         xlim([0 360])
%         ylim([-90 90])
%         title('log10 mean catch All L (g m^-^2)')
%         colormap('jet')
%         colorbar('h')
%         caxis([-5 -2])
%         stamp(cfile)
%         print('-dpng',[ppath 'Climatol_' harv '_global_AllL_catch.png'])
%
%% FracPFvD
% figure(27)
% surf(lon,lat,FracPFvD); view(2); hold on;
% shading flat
% title('(P+F):D mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPFvD.png'])
%
% % FracPDs
% figure(28)
% surf(lon,lat,FracPDs); view(2); hold on;
% shading flat
% title('SP:SD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPDs.png'])
%
% % FracPFs
% figure(29)
% surf(lon,lat,FracPFs); view(2); hold on;
% shading flat
% title('SP:SF mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPFs.png'])
%
% % FracPFvDs
% figure(30)
% surf(lon,lat,FracPFvDs); view(2); hold on;
% shading flat
% title('(SP+SF):SD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPFvDs.png'])
%
% % FracPDm
% figure(31)
% surf(lon,lat,FracPDm); view(2); hold on;
% shading flat
% title('MP:MD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPDm.png'])
%
% % FracPFm
% figure(32)
% surf(lon,lat,FracPFm); view(2); hold on;
% shading flat
% title('MP:MF mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPFm.png'])
%
% % FracPFvDm
% figure(33)
% surf(lon,lat,FracPFvDm); view(2); hold on;
% shading flat
% title('(MP+MF):MD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Climatol_' harv '_global_FracPFvDm.png'])
%


