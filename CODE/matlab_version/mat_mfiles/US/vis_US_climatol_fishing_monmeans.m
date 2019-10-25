% Visualize output of FEISTY in US LMEs
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
%cfile = 'Dc_enc50-b210_m4-b175-k060_c50-b250_D075_J075_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
%harv = 'fish_F030_P060_D060';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/US/'];
if (~isdir(ppath))
    mkdir(ppath)
end

%load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'Means_Climatol_' harv '_' cfile '.mat']);
%load([fpath 'Climatol_All_fish03_ICbiom_1e-10.mat']);

close all

%% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=10; %Set these bounds for your data
plotmaxlat=75;
plotminlon=-180;
plotmaxlon=-50;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;


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


%% bent
figure(50)
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(Zb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean Benthic inverts (g m^-^2)')
set(gca,'XColor','none','YColor','none')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_US_BENT.png'])

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
FracLM = AllL ./ (AllM+AllL);
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
% axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,log10(All))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 2]);
% hcb = colorbar('h');
% ylim(hcb,[-1 2])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All fishes (g m^-^2)')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_US_All.png'])
% 
% % all F
% figure(22)
% axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,log10(AllF))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% hcb = colorbar('h');
% ylim(hcb,[-1 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All F (g m^-^2)')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_US_AllF.png'])
% 
% % all D
% figure(23)
% axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,log10(AllD))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% hcb = colorbar('h');
% ylim(hcb,[-1 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All D (g m^-^2)')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_US_AllD.png'])
% 
% % All P
% figure(24)
% axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(geolat_t,geolon_t,log10(AllP))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-1 1]);
% hcb = colorbar('h');
% ylim(hcb,[-1 1])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Climatology log10 mean All P (g m^-^2)')
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Climatol_' harv '_US_AllP.png'])


%% All 4 on subplots
figure(27)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
%     hcb = colorbar('h');
set(gcf,'renderer','painters')
%title('log10 mean All F (g m^-^2)')
text(-1.3,1.85,'A','FontSize',14)
set(gca,'XColor','none','YColor','none')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
%     hcb = colorbar('h');
set(gcf,'renderer','painters')
%title('log10 mean All P (g m^-^2)')
text(-1.3,1.85,'B','FontSize',14)
set(gca,'XColor','none','YColor','none')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
%     hcb = colorbar('h');
set(gcf,'renderer','painters')
%title('log10 mean All D (g m^-^2)')
text(-1.3,1.85,'C','FontSize',14)
set(gca,'XColor','none','YColor','none')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
%     hcb = colorbar('h');
set(gcf,'renderer','painters')
colorbar('Position',[0.25 0.5 0.5 0.0255],'orientation','horizontal')
%title('log10 mean All fishes (g m^-^2)')
text(-1.3,1.85,'D','FontSize',14)
%     stamp([harv '_' cfile])
set(gca,'XColor','none','YColor','none')
print('-dpng',[ppath 'Climatol_' harv '_US_All_subplot.png'])


%% Ratios on subplots red-white-blue
figure(29)
% all P:F
subplot('Position',[0 0.55 1 0.4])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
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
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FracPD)
%     colormap(cmap_color_rb)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('B. Large Pelagics : Demersals')
%     stamp([harv '_' cfile])
%print('-dpng',[ppath 'Climatol_' harv '_US_ratios_subplot_v2.png'])

%% 3 figure subplot P:D, P:F, M:L
figure(30)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FracPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.2 0.485 0.6 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_US_ratios_subplot_v3.png'])


