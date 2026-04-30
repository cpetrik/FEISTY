% Plot differences between FEISTY only (2D) and COBALT-FEISTY (3D)
% Feisty biomasses

clear
close all

%%
gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
load([gpath 'Data_grid_OM4_05_COBALTv3.mat'],'GRD');
load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')
[ni,nj]=size(geolon);

dz = diff(z_l);

%% ONLINE -----------------------------------------------------------
npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

load([npath '19900101.ocean_feisty_tracers_z_means.mat'])

NtF = tSF+tMF;
NtP = tSP+tMP+tLP;
NtD = tSD+tMD+tLD;
NtB = tBE;
NtAll = NtF + NtP + NtD;

NsB = sBE;
NsF = sSF+sMF;
NsP = sSP+sMP+sLP;
NsD = sSD+sMD+sLD;
NsM = sMF+sMP+sMD;
NsL = sLP+sLD;
NsAll = NsF+NsP+NsD;
NFracPD = NsP ./ (NsP + NsD);
NFracPF = NsP ./ (NsP + NsF);
NFracLM = NsL ./ (NsL + NsM);

%% OFFLINE -----------------------------------------------------------
fpath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/COBALTv3_Hindcast_HalfDeg/';

vers = 'Historic_1990_All_fish03';
mod = 'Historic1990_All_fish03_2D_offline';

load([fpath 'Time_Means_' vers '_' cfile '.mat']);
load([fpath 'Annual_Means_' vers '_' cfile '.mat']);


FtSF = sf_tmean(1:12);
FtMF = mf_tmean(1:12);
FtSP = sp_tmean(1:12);
FtMP = mp_tmean(1:12);
FtLP = lp_tmean(1:12);
FtSD = sd_tmean(1:12);
FtMD = md_tmean(1:12);
FtLD = ld_tmean(1:12);

FtF = sf_tmean(1:12)+mf_tmean(1:12);
FtP = sp_tmean(1:12)+mp_tmean(1:12)+lp_tmean(1:12);
FtD = sd_tmean(1:12)+md_tmean(1:12)+ld_tmean(1:12);
FtB = b_tmean(1:12);
FtAll = FtF + FtP + FtD;

FZsf=NaN*ones(ni,nj);
FZsp=NaN*ones(ni,nj);
FZsd=NaN*ones(ni,nj);
FZmf=NaN*ones(ni,nj);
FZmp=NaN*ones(ni,nj);
FZmd=NaN*ones(ni,nj);
FZlp=NaN*ones(ni,nj);
FZld=NaN*ones(ni,nj);
FZb=NaN*ones(ni,nj);

FZsf(GRD.ID)=sf_abio(:,1);
FZsp(GRD.ID)=sp_abio(:,1);
FZsd(GRD.ID)=sd_abio(:,1);
FZmf(GRD.ID)=mf_abio(:,1);
FZmp(GRD.ID)=mp_abio(:,1);
FZmd(GRD.ID)=md_abio(:,1);
FZlp(GRD.ID)=lp_abio(:,1);
FZld(GRD.ID)=ld_abio(:,1);
FZb(GRD.ID)=b_abio(:,1);

% Diff maps of all fish
FsAll = FZsp+FZsf+FZsd+FZmp+FZmf+FZmd+FZlp+FZld;
FsF = FZsf+FZmf;
FsP = FZsp+FZmp+FZlp;
FsD = FZsd+FZmd+FZld;
FsM = FZmp+FZmf+FZmd;
FsL = FZlp+FZld;
FFracPD = FsP ./ (FsP + FsD);
FFracPF = FsP ./ (FsP + FsF);
FFracLM = FsL ./ (FsL + FsM);


%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    

set(groot,'defaultAxesColorOrder',cm10);

%%
geolon = double(geolon);
geolat = double(geolat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% Diffs Online - Offline

dtB = NtB - FtB;
dtF = NtF - FtF;
dtP = NtP - FtP;
dtD = NtD - FtD;
dtAll = NtAll - FtAll;

dsB = NsB - FZb;
dsF = NsF - FsF;
dsP = NsP - FsP;
dsD = NsD - FsD;
dsAll = NsAll - FsAll;

%% PLOTS TOGETHER
% Time series
tmos = 1:12;

% Log10
figure(1)
plot(tmos,log10(NtF(tmos)+eps),'r','LineWidth',2); hold on;
plot(tmos,log10(NtP(tmos)+eps),'b','LineWidth',2); hold on;
plot(tmos,log10(NtD(tmos)+eps),'color',cm10(2,:),'LineWidth',2); hold on;
plot(tmos,log10(FtF(tmos)+eps),'--r','LineWidth',2); hold on;
plot(tmos,log10(FtP(tmos)+eps),'--b','LineWidth',2); hold on;
plot(tmos,log10(FtD(tmos)+eps),'--','color',cm10(2,:),'LineWidth',2); hold on;
legend({'OnF','OnP','OnD','OffF','OffP','OffD'})
legend('location','west')
title('log_1_0 Integrated Biomass (gC m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_ts_log10mean_feisty_types.png'])

figure(2)
subplot(2,2,1)
plot(tmos,(NtF(tmos)+eps),'r','LineWidth',2); hold on;
plot(tmos,(FtF(tmos)+eps),'--r','LineWidth',2); hold on;
title('Forage')
subplot(2,2,2)
plot(tmos,(NtP(tmos)+eps),'b','LineWidth',2); hold on;
plot(tmos,(FtP(tmos)+eps),'--b','LineWidth',2); hold on;
title('Lg Pelagics')
subplot(2,2,3)
plot(tmos,(NtD(tmos)+eps),'color',cm10(2,:),'LineWidth',2); hold on;
plot(tmos,(FtD(tmos)+eps),'--','color',cm10(2,:),'LineWidth',2); hold on;
title('Demersals')
subplot(2,2,4)
plot(tmos,(NtAll(tmos)+eps),'k','LineWidth',2); hold on;
plot(tmos,(FtAll(tmos)+eps),'--k','LineWidth',2); hold on;
title('All fishes (gC m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_ts_mean_feisty_types.png'])

%% Maps
% Forage
figure(3)
subplot(2,1,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsF)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(0.5,-4,'log10 mean Forage (gC m^-^2)','HorizontalAlignment','center')

subplot(2,1,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsF)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_map_F.png'])

%%
figure(4)
subplot(2,1,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsP)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(0.5,-4,'log10 mean Lg Pelagic (gC m^-^2)','HorizontalAlignment','center')

subplot(2,1,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsP)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_map_P.png'])

% Dem
figure(5)
subplot(2,1,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsD)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(0.5,-4,'log10 mean Demersal (gC m^-^2)','HorizontalAlignment','center')

subplot(2,1,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsD)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_map_D.png'])

%% All
figure(6)
subplot(2,1,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsAll)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(0.5,-4,'log10 mean All fishes (gC m^-^2)','HorizontalAlignment','center')

subplot(2,1,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsAll)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_map_AllFishes.png'])

% Bent
figure(7)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FZb)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(0.5,-4,'log10 mean Benthos (gC m^-^2)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsB)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_map_mean_Bent.png'])


%% PLOTS OF DIFFS
% Time series
% NO Log10
figure(8)
plot(tmos,(dtF(tmos)+eps),'r','LineWidth',2); hold on;
plot(tmos,(dtP(tmos)+eps),'b','LineWidth',2); hold on;
plot(tmos,(dtD(tmos)+eps),'color',cm10(2,:),'LineWidth',2); hold on; 
legend({'F','P','D'})
legend('location','west')
title('Diff Integrated Biomass (gC m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_ts_diff_feisty_types.png'])

%% Maps
LAT = geolat;
LON= geolon;
figure(9)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,dsF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-20 20]);
colorbar
set(gcf,'renderer','painters')
title('Diff mean All F (g m^-^2)')
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,dsD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
colorbar
set(gcf,'renderer','painters')
title('Diff mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,dsP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-150 150]);
colorbar
set(gcf,'renderer','painters')
title('Diff mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,dsAll)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-150 150]);
colorbar
set(gcf,'renderer','painters')
title('Diff mean All fishes (g m^-^2)')
stamp('')
print('-dpng',[ppath mod '_OM4_05_COBALTv3_FEISTY_on3D-off2D_global_diffs.png'])

%% Bent
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsB)))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Benthic biomass difference (gC m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on3D-off2D_map_mean_diff_Bent.png'])

