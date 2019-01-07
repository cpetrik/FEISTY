% Visualize output of FEISTY in Gulf of Mexico
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/GoMex/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);
load([fpath 'GoMex_ts_means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

close all

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %15; %Set these bounds for your data
plotmaxlat=90; %35;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[100W 75W] = GoMex

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


%% Plots in time
y = time;
nt = length(time);

% All size classes of all

figure(4)
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title(['Climatol ' tharv])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_all_sizes.png'])

figure(5)
F = sf_tmean(1:nt)+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;

plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['Climatol ' tharv])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_all_types.png'])

% FISHING All size classes of all

figure(6)
plot(y,log10(mf_tmy),'color',[0 0.7 0],'Linewidth',1); hold on;
plot(y,log10(mp_tmy),'color',[1 0 0],'Linewidth',1); hold on;
plot(y,log10(lp_tmy),'color',[0.5 0 0],'Linewidth',1); hold on;
plot(y,log10(md_tmy),'color',[0 0.5 0.75],'Linewidth',1); hold on;
plot(y,log10(ld_tmy),'color',[0 0 0.75],'Linewidth',1); hold on;
legend('MF','MP','LP','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-7 0])
xlabel('Time (mo)')
ylabel('log10 Catch (g m^-^2)')
title(['Climatol ' tharv])
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_catch_all_sizes.png'])

figure(7)
F = mf_tmy;
P = mp_tmy+lp_tmy;
D = md_tmy+ld_tmy;

plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-7 0])
xlabel('Time (y)')
ylabel('log10 Catch (g m^-^2)')
title(['Climatol ' tharv])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_catch_all_types.png'])


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

% ocean cells
% figure(55)
% surf(lon,lat,ocean); view(2); hold on;
% shading flat
% title('Water cells')
% colormap('jet')
% colorbar('h')
% caxis([1 2])
% stamp([harv '_' cfile])
% print('-dpng',[ppath 'Ocean_cells.png'])

%%
latlim=[15 35];
lonlim=[-100 -75];

% bent
figure(50)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
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
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_BENT.png'])

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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_SP.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_SF.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_SD.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_MP.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_MF.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_MD.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_LP.png'])
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
% print('-dpng',[ppath 'Climatol_' harv '_GoMex_LD.png'])

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

%% ALL
figure(21)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
hcb = colorbar('h');
ylim(hcb,[-1 2])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 mean All fishes (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_All.png'])

% all F
figure(22)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 mean All F (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_AllF.png'])

% all D
figure(23)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 mean All D (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_AllD.png'])

% All P
figure(24)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology log10 mean All P (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_AllP.png'])

%% All 4 on subplots
figure(27)
% all F
subplot('Position',[0 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Forage fishes (log_1_0 g m^-^2)')

% all D
subplot('Position',[0 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
%     hcb = colorbar('h');
%     ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Demersals (log_1_0 g m^-^2)')

% All P
subplot('Position',[0.48 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('Large pelagics (log_1_0 g m^-^2)')

% All
subplot('Position',[0.48 0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.9 0.25 0.03 0.5],'orientation','vertical')
set(gcf,'renderer','painters')
%title('log10 mean All fishes (g m^-^2)')
title('All fishes (log_1_0 g m^-^2)')
%     stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_All_subplot.png'])

%% 3 figure subplot P:D, P:F, M:L
figure(30)
subplot('Position',[0 0.5 0.48 0.45])
%P:D
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.5 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPF)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.48 0.45])
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracLM)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.7 0.02 0.035 0.425],'orientation','vertical')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_ratios_subplot_v3.png'])


%% Interpolate to a finer grid
%geolon_t: 0.5:359.5
%geolat_t: -89.5:89.5
lats = -89.5:0.25:89.5;
lons = 0.5:0.25:359.5;
[glon,glat] = meshgrid(lons,lats);


hF = griddata(geolat_t,geolon_t,AllF,glat,glon);
hP = griddata(geolat_t,geolon_t,AllP,glat,glon);
hD = griddata(geolat_t,geolon_t,AllD,glat,glon);
hS = griddata(geolat_t,geolon_t,AllS,glat,glon);
hM = griddata(geolat_t,geolon_t,AllM,glat,glon);
hL = griddata(geolat_t,geolon_t,AllL,glat,glon);

%%
figure(100)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean F (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_F_highres.png'])

figure(101)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean P (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_P_highres.png'])

figure(102)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Climatology log10 mean D (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_GoMex_D_highres.png'])



