% nuts, LTL, and feisty forcing
% means by latitude

clear
close all

%% molN/kg to gC/m3 or gC/m3/d
S2D= 60*60*24;
D2Y = S2D*365;
%MZ and LZ molN/kg
N2Ckg= 1035 * (106/16) * 12.01;
%HPloss 'mol N kg-1 s-1'
N2CkgD = N2Ckg * S2D;
N2CkgY = N2Ckg * D2Y;
%fn_tot_btm = 'mol m-2 s-1'
N2CmD = (106/16) * 12.01 * S2D;
N2CmY = (106/16) * 12.01 * D2Y;

% molN/kg to gC/m3
NtoC= 1035 * (106/16) * 12.01;

%% ONLINE -----------------------------------------------------------
npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

load([npath 'ocean_cobalt_feisty_forcing_z.199001-199412_means_lat.mat'])

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
NmBDE = N2CmD * mDE;
NmMZ = N2Ckg * mMZ;
NmLZ = N2Ckg * mLZ;

NvMZ = N2Ckg * vMZ;
NvLZ = N2Ckg * vLZ;

clear mDE mMZ mLZ 
clear vMZ vLZ 

%% OFFLINE -----------------------------------------------------------
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412_means_lat.mat'])

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
FmBDE = N2CmD * mDE;
FmMZ = N2Ckg * mMZ;
FmLZ = N2Ckg * mLZ;

FvMZ = N2Ckg * vMZ;
FvLZ = N2Ckg * vLZ;

clear mDE mMZ mLZ 
clear vMZ vLZ 

%% COBALT PLANKTON TRACERS ---------------------------------------

%% ONLINE -----------------------------------------------------------
load([npath 'ocean_cobalt_tracers_month_z.199001-199412_means_lat.mat'])

%% Put everything in gC m^-^2 
NtBDE = NtoC * tDE(1:60);
NtDI = NtoC * tDI(1:60);
NtSP = NtoC * tSP(1:60);
NtLP = NtoC * tLP(1:60);
NtSZ = NtoC * tSZ(1:60);

NvDE = NtoC * vDE(:,1);
NvDI = NtoC * vDI(:,1);
NvSP = NtoC * vSP(:,1);
NvLP = NtoC * vLP(:,1);
NvSZ = NtoC * vSZ(:,1);

NsBDE = NtoC * sDE(:,:,1);
NsDI = NtoC * sDI(:,:,1);
NsSP = NtoC * sSP(:,:,1);
NsLP = NtoC * sLP(:,:,1);
NsSZ = NtoC * sSZ(:,:,1);

clear tDE tDI vDE vDI sDE sDI
clear tSP tLP tSZ vSP vLP vSZ sSP sLP sSZ

%% OFFLINE -----------------------------------------------------------
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_tracers_month_z.199001-199412_means_lat.mat'])

%% Put everything in gC m^-^2 
FtDE = NtoC * tDE(1:60);
FtDI = NtoC * tDI(1:60);
FtSP = NtoC * tSP(1:60);
FtLP = NtoC * tLP(1:60);
FtSZ = NtoC * tSZ(1:60);

FvDE = NtoC * vDE(:,1);
FvDI = NtoC * vDI(:,1);
FvSP = NtoC * vSP(:,1);
FvLP = NtoC * vLP(:,1);
FvSZ = NtoC * vSZ(:,1);

FsDE = NtoC * sDE(:,:,1);
FsDI = NtoC * sDI(:,:,1);
FsSP = NtoC * sSP(:,:,1);
FsLP = NtoC * sLP(:,:,1);
FsSZ = NtoC * sSZ(:,:,1);

clear tDE tDI vDE vDI sDE sDI
clear tSP tLP tSZ vSP vLP vSZ sSP sLP sSZ

%% COBALT NUT TRACERS ---------------------------------------

%%
load([fpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')

dz = diff(z_l);

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
[ni,nj]=size(geolon);
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

dtDE = NtBDE - FtDE;
dtDI = NtDI - FtDI;
dtSP = NtSP - FtSP;
dtLP = NtLP - FtLP;
dtSZ = NtSZ - FtSZ;

dsDE = NsBDE - FsDE;
dsDI = NsDI - FsDI;
dsSP = NsSP - FsSP;
dsLP = NsLP - FsLP;
dsSZ = NsSZ - FsSZ;

dvDE = NvDE - FvDE;
dvDI = NvDI - FvDI;
dvSP = NvSP - FvSP;
dvLP = NvLP - FvLP;
dvSZ = NvSZ - FvSZ;

%% PLOTS TOGETHER
% Time series
tmos = 1:60;

% Log10
figure(1)
subplot(3,1,1)
plot(tmos,log10(NtSP(tmos)+eps),'color',cm10(3,:)); hold on;
plot(tmos,log10(NtLP(tmos)+eps),'color',cm10(4,:)); hold on; 
plot(tmos,log10(NtDI(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,log10(FtSP(tmos)+eps),'--','color',cm10(3,:)); hold on;
plot(tmos,log10(FtLP(tmos)+eps),'--','color',cm10(4,:)); hold on; 
plot(tmos,log10(FtDI(tmos)+eps),'--','color',cm10(5,:)); hold on;
legend({'OnSP','OnLP','OnDI','OffSP','OffLP','OffDI'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
plot(tmos,log10(NtSZ(tmos)+eps),'color',cm10(8,:)); hold on;
plot(tmos,log10(FtSZ(tmos)+eps),'--','color',cm10(8,:)); hold on;
legend({'OnSZ','OffSZ'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^2)')
% xlabel('Months')
subplot(3,1,3)
plot(tmos,log10(NtBDE(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(FtDE(tmos)+eps),'--','color',cm10(1,:)); hold on;
legend({'OnDet','OffDet'})
legend('location','eastoutside')
title('log_1_0 Integrated Detritus (gC m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_log10mean_cobalt_tracers.png'])

%% No log
figure(2)
subplot(3,1,1)
plot(tmos,(NtSP(tmos)+eps),'color',cm10(3,:)); hold on;
plot(tmos,(NtLP(tmos)+eps),'color',cm10(4,:)); hold on; 
plot(tmos,(NtDI(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,(FtSP(tmos)+eps),'--','color',cm10(3,:)); hold on;
plot(tmos,(FtLP(tmos)+eps),'--','color',cm10(4,:)); hold on; 
plot(tmos,(FtDI(tmos)+eps),'--','color',cm10(5,:)); hold on;
legend({'OnSP','OnLP','OnDI','OffSP','OffLP','OffDI'})
legend('location','eastoutside')
title('Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
plot(tmos,(NtSZ(tmos)+eps),'color',cm10(8,:)); hold on;
plot(tmos,(FtSZ(tmos)+eps),'--','color',cm10(8,:)); hold on;
legend({'OnSZ','OffSZ'})
legend('location','eastoutside')
title('Integrated Biomass (gC m^-^2)')
% xlabel('Months')
subplot(3,1,3)
plot(tmos,(NtBDE(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,(FtDE(tmos)+eps),'--','color',cm10(1,:)); hold on;
legend({'OnDet','OffDet'})
legend('location','eastoutside')
title('Integrated Detritus (gC m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_mean_cobalt_tracers.png'])

%% Vert
figure(3)
subplot(2,2,3)
plot(log10(NvSP(:,1)+eps),-1*z_l,'color',cm10(3,:)); hold on;
plot(log10(NvLP(:,1)+eps),-1*z_l,'color',cm10(4,:)); hold on; 
plot(log10(NvDI(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on;
plot(log10(FvSP(:,1)+eps),-1*z_l,'--','color',cm10(3,:)); hold on;
plot(log10(FvLP(:,1)+eps),-1*z_l,'--','color',cm10(4,:)); hold on; 
plot(log10(FvDI(:,1)+eps),-1*z_l,'--','color',cm10(5,:)); hold on;
title('log_1_0')
subplot(2,2,4)
plot(log10(NvSZ(:,1)+eps),-1*z_l,'color',cm10(8,:)); hold on;
plot(log10(FvSZ(:,1)+eps),-1*z_l,'--','color',cm10(8,:)); hold on;
title('log_1_0')
ylabel('Depth (m)')

subplot(2,2,1)
plot((NvSP(1:8,1)),-1*z_l(1:8),'color',cm10(3,:)); hold on;
plot((NvLP(1:8,1)),-1*z_l(1:8),'color',cm10(4,:)); hold on; 
plot((NvDI(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on;
plot((FvSP(1:8,1)),-1*z_l(1:8),'--','color',cm10(3,:)); hold on;
plot((FvLP(1:8,1)),-1*z_l(1:8),'--','color',cm10(4,:)); hold on; 
plot((FvDI(1:8,1)),-1*z_l(1:8),'--','color',cm10(5,:)); hold on;
legend({'OnSP','OnLP','OnDI','OffSP','OffLP','OffDI'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
subplot(2,2,2)
plot((NvSZ(1:8,1)),-1*z_l(1:8),'color',cm10(8,:)); hold on;
plot((FvSZ(1:8,1)),-1*z_l(1:8),'--','color',cm10(8,:)); hold on;
legend({'OnSZ','OffSZ'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_mean_cobalt_tracers.png'])

%% Maps
% Phyto bio 
figure(4)
subplot(2,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsSP(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean Small phyto (gC m^-^2)','HorizontalAlignment','center')

subplot(2,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsSP(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')

subplot(2,2,3) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsLP(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean Large phyto (gC m^-^2)','HorizontalAlignment','center')

subplot(2,2,4) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsLP(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_sm_lg_phyto.png'])


%% SZ
figure(5)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsSZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.5 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean Small zoo (gC m^-^2)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsSZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.5 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_SZ.png'])

%% DIAZ
figure(6)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsDI(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean Diaz (gC m^-^2)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsDI(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Diaz.png'])


%% Det 
figure(17)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsDE(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean btm detritus flux (gC m^-^2 d^-^1)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsBDE(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_Det.png'])


%% PLOTS OF DIFFS ---------------------------------------
% Time series
% NO Log10
figure(7)
subplot(3,1,1)
plot(tmos,(dtSP(tmos)+eps),'color',cm10(3,:),'LineWidth',2); hold on;
plot(tmos,(dtLP(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on; 
plot(tmos,(dtDI(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
title('Integrated Biomass (gC m^-^2) Difference Online - Offline')
legend({'SP','LP','Diaz'})
legend('location','southwest')
subplot(3,1,2)
plot(tmos,(dtSZ(tmos)+eps),'color',cm10(8,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated SZ Biomass (gC m^-^2) Difference Online - Offline')
% xlabel('Months')
subplot(3,1,3)
plot(tmos,(dtDE(tmos)+eps),'color',cm10(1,:),'LineWidth',2); 
title('Detritus Concentration (gC m^-^2) Difference Online - Offline')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_diff_cobalt_tracers.png'])

%% Vert
figure(8)
subplot(1,2,1)
plot((dvSP(1:9,1)),-1*z_l(1:9),'color',cm10(3,:),'LineWidth',2); hold on;
plot((dvLP(1:9,1)),-1*z_l(1:9),'color',cm10(4,:),'LineWidth',2); hold on; 
plot((dvDI(1:9,1)),-1*z_l(1:9),'color',cm10(5,:),'LineWidth',2); hold on;
legend({'SP','LP','Diaz'})
legend('location','southwest')
title('Mean Biomass (gC m^-^3) Difference Online - Offline')
subplot(1,2,2)
plot((dvSZ(1:9,1)),-1*z_l(1:9),'color',cm10(8,:),'LineWidth',2); hold on;
title('Mean Biomass (gC m^-^3) Difference Online - Offline')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_diff_cobalt_tracers.png'])


%% Maps
% DI, SP, LP
figure(9)
subplot(2,2,1) %SP
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsSP(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Small phyto biomass difference (gC m^2)')
text(6,2.5,'Online - Offline','HorizontalAlignment','center')

subplot(2,2,2) %LP
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLP(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Large phyto biomass difference (gC m^2)')

subplot(2,2,3) %DI
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDI(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Diazotroph biomass difference (gC m^2)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_phytos.png'])

%% SZ
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsSZ(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Small zoo difference (gC m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_SZ.png'])

%% Det
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDE(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Detrius difference (gC m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_Det.png'])

%% SP
figure(12) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsSP(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Small phyto biomass difference (gC m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_SP.png'])

%% LP
figure(13)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLP(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Large phyto biomass difference (gC m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_LP.png'])

%% DI
figure(14)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDI(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Diazotroph biomass difference (gC m^2)')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_Diaz.png'])







%% COBALT NUTS  ---------------------------------------




%% molN/kg to gC/m3
NtoC= 1035 * (106/16) * 12.01;

%% ONLINE -----------------------------------------------------------
%npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';
npath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FEISTYon/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

load([npath 'ocean_cobalt_nuts_month_z.199001-199412_means.mat'])

%% Nothing needs to be in gC m^-^2 
% Should it be grams N?
% Chl from ug/kg to ug/m3
NtNH4 = tNH4(1:60);
NtNO3 = tNO3(1:60);
NtO2  = tO2(1:60);
NtCHL = tCHL(1:60)*1035;

NvNH4 = vNH4(:,1);
NvNO3 = vNO3(:,1);
NvO2  = vO2(:,1);
NvCHL = vCHL(:,1)*1035;

NsNH4 = sNH4(:,:,1);
NsNO3 = sNO3(:,:,1);
NsO2 = sO2(:,:,1);
NsCHL = sCHL(:,:,1)*1035;
NsCHLs = sCHLs(:,:,1)*1035;

clear tNH4 tNO3 vNH4 vNO3 sNH4 sNO3
clear tO2 tCHL vO2 vCHL sO2 sCHL sCHLs

%% OFFLINE -----------------------------------------------------------
%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
fpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/COBALTonly/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_nuts_month_z.199001-199412_means.mat'])

%% Put everything in gC m^-^2 
FtNH4 = tNH4(1:60);
FtNO3 = tNO3(1:60);
FtO2  = tO2(1:60);
FtCHL = tCHL(1:60)*1035;

FvNH4 = vNH4(:,1);
FvNO3 = vNO3(:,1);
FvO2  = vO2(:,1);
FvCHL = vCHL(:,1)*1035;

FsNH4 = sNH4(:,:,1);
FsNO3 = sNO3(:,:,1);
FsO2  = sO2(:,:,1);
FsCHL = sCHL(:,:,1)*1035;
FsCHLs = sCHLs(:,:,1)*1035;

clear tNH4 tNO3 vNH4 vNO3 sNH4 sNO3
clear tO2 tCHL vO2 vCHL sO2 sCHL sCHLs

FsNO3(FsNO3<0) = 0;
FsO2(FsO2<0) = 0;
NsO2(NsO2<0) = 0;

%%
load([fpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')

dz = diff(z_l);

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
[ni,nj]=size(geolon);
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

dtNH4 = NtNH4 - FtNH4;
dtNO3 = NtNO3 - FtNO3;
dtO2 = NtO2 - FtO2;
dtCHL = NtCHL - FtCHL;

dsNH4 = NsNH4 - FsNH4;
dsNO3 = NsNO3 - FsNO3;
dsO2 = NsO2 - FsO2;
dsCHL = NsCHL - FsCHL;
dsCHLs = NsCHLs - FsCHLs;

dvNH4 = NvNH4 - FvNH4;
dvNO3 = NvNO3 - FvNO3;
dvO2 = NvO2 - FvO2;
dvCHL = NvCHL - FvCHL;

%% PLOTS TOGETHER
% Time series
tmos = 1:60;

% Log10
figure(1)
subplot(4,1,1)
plot(tmos,log10(NtNO3(tmos)+eps),'color',cm10(10,:)); hold on;
plot(tmos,log10(FtNO3(tmos)+eps),'--','color',cm10(10,:)); hold on;
legend({'OnNO3','OffNO3'})
legend('location','eastoutside')
title('log_1_0 Integrated Concentration (mol m^-^2)')
subplot(4,1,2)
plot(tmos,log10(NtNH4(tmos)+eps),'color',cm10(9,:)); hold on;
plot(tmos,log10(FtNH4(tmos)+eps),'--','color',cm10(9,:)); hold on;
legend({'OnNH4','OffNH4'})
legend('location','eastoutside')
title('log_1_0 Integrated Concentration (mol m^-^2)')
subplot(4,1,3)
plot(tmos,log10(1e-3*NtCHL(tmos)+eps),'color',cm10(2,:)); hold on;
plot(tmos,log10(1e-3*FtCHL(tmos)+eps),'--','color',cm10(2,:)); hold on;
legend({'OnChl','OffChl'})
legend('location','eastoutside')
title('log_1_0 Integrated Chl (mg m^-^2)')
%ylim([4.45 4.75])
subplot(4,1,4)
plot(tmos,log10(NtO2(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(FtO2(tmos)+eps),'--','color',cm10(1,:)); hold on;
legend({'OnO2','OffO2'})
legend('location','eastoutside')
title('log_1_0 Integrated O2 (mol m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_log10mean_cobalt_nuts.png'])

%% No log
figure(2)
subplot(4,1,1)
plot(tmos,(NtNO3(tmos)+eps),'color',cm10(10,:)); hold on;
plot(tmos,(FtNO3(tmos)+eps),'--','color',cm10(10,:)); hold on;
legend({'OnNO3','OffNO3'})
legend('location','eastoutside')
title('Integrated Concentration (mol m^-^2)')
subplot(4,1,2)
plot(tmos,(NtNH4(tmos)+eps),'color',cm10(9,:)); hold on;
plot(tmos,(FtNH4(tmos)+eps),'--','color',cm10(9,:)); hold on;
legend({'OnNH4','OffNH4'})
legend('location','eastoutside')
title('Integrated Concentration (mol m^-^2)')
subplot(4,1,3)
plot(tmos,(1e-3*NtCHL(tmos)+eps),'color',cm10(2,:)); hold on;
plot(tmos,(1e-3*FtCHL(tmos)+eps),'--','color',cm10(2,:)); hold on;
legend({'OnChl','OffChl'})
legend('location','eastoutside')
title('Integrated Chl (mg m^-^2)')
subplot(4,1,4)
plot(tmos,(NtO2(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,(FtO2(tmos)+eps),'--','color',cm10(1,:)); hold on;
legend({'OnO2','OffO2'})
legend('location','eastoutside')
title('Integrated O2 (mol m^-^2)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_mean_cobalt_nuts.png'])

%% Vert
figure(3)
subplot(2,3,4)
plot(log10(NvNH4(:,1)+eps),-1*z_l,'color',cm10(9,:)); hold on;
plot(log10(NvNO3(:,1)+eps),-1*z_l,'color',cm10(10,:)); hold on;
plot(log10(FvNH4(:,1)+eps),-1*z_l,'--','color',cm10(9,:)); hold on;
plot(log10(FvNO3(:,1)+eps),-1*z_l,'--','color',cm10(10,:)); hold on;
title('log_1_0 N')
subplot(2,3,5)
plot(log10(NvO2(:,1)+eps),-1*z_l,'color',cm10(1,:)); hold on;
plot(log10(FvO2(:,1)+eps),-1*z_l,'--','color',cm10(1,:)); hold on;
title('log_1_0 O2')
subplot(2,3,6)
plot(log10(NvCHL(:,1)+eps),-1*z_l,'color',cm10(2,:)); hold on;
plot(log10(FvCHL(:,1)+eps),-1*z_l,'--','color',cm10(2,:)); hold on;
title('log_1_0 Chl')
ylabel('Depth (m)')

subplot(2,3,1)
plot((NvNH4(1:8,1)),-1*z_l(1:8),'color',cm10(9,:)); hold on; 
plot((NvNO3(1:8,1)),-1*z_l(1:8),'color',cm10(10,:)); hold on;
plot((FvNH4(1:8,1)),-1*z_l(1:8),'--','color',cm10(9,:)); hold on; 
plot((FvNO3(1:8,1)),-1*z_l(1:8),'--','color',cm10(10,:)); hold on;
legend({'OnNH4','OnNO3','OffNH4','OffNO3'})
legend('location','southwest')
title('Mean N concentration (mol m^-^3)')
subplot(2,3,2)
plot((NvO2(1:8,1)),-1*z_l(1:8),'color',cm10(1,:)); hold on;
plot((FvO2(1:8,1)),-1*z_l(1:8),'--','color',cm10(1,:)); hold on;
legend({'OnO2','OffO2'})
legend('location','southeast')
title('Mean O2 concentration (mol m^-^3)')
subplot(2,3,3)
plot((NvCHL(1:8,1)),-1*z_l(1:8),'color',cm10(2,:)); hold on;
plot((FvCHL(1:8,1)),-1*z_l(1:8),'--','color',cm10(2,:)); hold on;
legend({'OnChl','OffChl'})
legend('location','southeast')
title('Mean Chl concentraion (ug m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_mean_cobalt_nuts.png'])

%% Maps
% Ns
figure(4)
subplot(2,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsNO3(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -4]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean NO3 (mol m^-^2)','HorizontalAlignment','center')

subplot(2,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsNO3(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -4]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')

subplot(2,2,3) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsNH4(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-8 -6]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean NH4 (mol m^-^2)','HorizontalAlignment','center')

subplot(2,2,4) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsNH4(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-8 -6]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_no3_nh4.png'])


%% CHLs
figure(5)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsCHLs(:,:,1))))
cmocean('algae')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([1 4]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean surf Chl (ug m^-^3)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsCHLs(:,:,1))))
cmocean('algae')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([1 4]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_surfChl.png'])

%% Chl int
figure(6)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(1e-3*FsCHL(:,:,1))))
cmocean('algae')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([1 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean integrated Chl (mg m^-^2)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(1e-3*NsCHL(:,:,1))))
cmocean('algae')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([1 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_intChl.png'])


%% O2
figure(17)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsO2(:,:,1))))
cmocean('grey')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor','w'); %[0.75 0.75 0.75]);
clim([-4.5 -3.25]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean O2 (mol m^-^2)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsO2(:,:,1))))
cmocean('grey')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor','w'); %,[0.75 0.75 0.75]);
clim([-4.5 -3.25]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_O2.png'])


%% PLOTS OF DIFFS ---------------------------------------
% Time series
% NO Log10
figure(7)
subplot(3,1,1)
plot(tmos,(dtNH4(tmos)+eps),'color',cm10(9,:),'LineWidth',2); hold on;
plot(tmos,(dtNO3(tmos)+eps),'color',cm10(10,:),'LineWidth',2); hold on;
title('Integrated Concentration (mol m^-^2) Difference Online - Offline')
legend({'NH4','NO3'})
legend('location','northwest')
subplot(3,1,2)
plot(tmos,(dtCHL(tmos)+eps),'color',cm10(2,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated Chl Concentration (m^-^2) Difference Online - Offline')
subplot(3,1,3)
plot(tmos,(dtO2(tmos)+eps),'color',cm10(1,:),'LineWidth',2); 
title('O2 Concentration (mol m^-^2) Difference Online - Offline')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_diff_cobalt_nuts.png'])

%% Vert
%CHl
figure(8)
subplot(1,2,1)
plot((dvNH4(1:9,1)),-1*z_l(1:9),'color',cm10(9,:),'LineWidth',2); hold on; 
plot((dvNO3(1:9,1)),-1*z_l(1:9),'color',cm10(10,:),'LineWidth',2); hold on;
plot((dvO2(1:9,1)),-1*z_l(1:9),'color',cm10(1,:),'LineWidth',2); hold on;
legend({'NH4','NO3','O2'})
legend('location','northeast')
title('Mean Concen (mol m^-^3) Difference Online - Offline')
subplot(1,2,2)
plot((dvCHL(1:9,1)),-1*z_l(1:9),'color',cm10(2,:),'LineWidth',2); hold on;
title('Mean Chl Concen (ug m^-^3) Difference Online - Offline')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_diff_cobalt_nuts.png'])


%% Maps
% NO3, NH4
figure(9)
subplot(1,2,1) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNO3(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3e-6 3e-6]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NO3 difference (mol m^2)')
text(4,2.5,'Online - Offline','HorizontalAlignment','center')

subplot(1,2,2) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNH4(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1e-7 1e-7]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NH4 difference (mol m^2)')

print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_no3_nh4.png'])

%% CHLs
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(1e-3*dsCHLs(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Surf Chl difference (mg m^3) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_surfChl.png'])

%% Chl
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(1e-3*dsCHL(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Int Chl difference (mg m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_intChl.png'])

%% O2
figure(12) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsO2(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5e-6 5e-6]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('O2 difference (mol m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_O2.png'])

%% NH4
figure(13)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNH4(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1e-7 1e-7]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NH4 difference (mol m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_NH4.png'])

%% NO3
figure(14)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNO3(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3e-6 3e-6]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NO3 difference (mol m^2)')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_NO3.png'])




%%
load([fpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')

dz = diff(z_l);

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
[ni,nj]=size(geolon);
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

dtDE = NtBDE - FtDE;
dtMZ = NtMZ - FtMZ;
dtLZ = NtLZ - FtLZ;
dtMH = NtMH - FtMH;
dtLH = NtLH - FtLH;

dsDE = NsBDE - FsDE;
dsMZ = NsMZ - FsMZ;
dsLZ = NsLZ - FsLZ;
dsMH = NsMH - FsMH;
dsLH = NsLH - FsLH;

dvMZ = NvMZ - FvMZ;
dvLZ = NvLZ - FvLZ;
dvMH = NvMH - FvMH;
dvLH = NvLH - FvLH;

%% PLOTS TOGETHER
% Time series
tmos = 1:120;

% Log10
figure(1)
subplot(3,1,1)
plot(tmos,log10(NtMZ(tmos)+eps),'color',cm10(4,:)); hold on;
plot(tmos,log10(NtLZ(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,log10(FtMZ(tmos)+eps),'--','color',cm10(4,:)); hold on;
plot(tmos,log10(FtLZ(tmos)+eps),'--','color',cm10(5,:)); hold on;
legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
plot(tmos,log10(NtMH(tmos)+eps),'color',cm10(6,:)); hold on; 
plot(tmos,log10(NtLH(tmos)+eps),'color',cm10(7,:)); hold on;
plot(tmos,log10(FtMH(tmos)+eps),'--','color',cm10(6,:)); hold on; 
plot(tmos,log10(FtLH(tmos)+eps),'--','color',cm10(7,:)); hold on;
legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
legend('location','eastoutside')
title('log_1_0 Integrated Higher Predation Rate (gC m^-^2 d^-^1)')
subplot(3,1,3)
plot(tmos,log10(NtBDE(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(FtDE(tmos)+eps),'--','color',cm10(1,:)); hold on;
legend({'OnBtmDet','OffBtmDet'})
legend('location','eastoutside')
title('log_1_0 Bottom Detritus Flux (gC m^-^2 d^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_log10mean_feisty_forcing.png'])

% No log
figure(2)
subplot(3,1,1)
plot(tmos,(NtMZ(tmos)+eps),'color',cm10(4,:)); hold on;
plot(tmos,(NtLZ(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,(FtMZ(tmos)+eps),'--','color',cm10(4,:)); hold on;
plot(tmos,(FtLZ(tmos)+eps),'--','color',cm10(5,:)); hold on;
legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
legend('location','eastoutside')
title('Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
plot(tmos,(NtMH(tmos)+eps),'color',cm10(6,:)); hold on; 
plot(tmos,(NtLH(tmos)+eps),'color',cm10(7,:)); hold on;
plot(tmos,(FtMH(tmos)+eps),'--','color',cm10(6,:)); hold on; 
plot(tmos,(FtLH(tmos)+eps),'--','color',cm10(7,:)); hold on;
legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
legend('location','eastoutside')
title('Integrated Higher Predation Rate (gC m^-^2 d^-^1)')
subplot(3,1,3)
plot(tmos,(NtBDE(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,(FtDE(tmos)+eps),'--','color',cm10(1,:)); hold on;
legend({'OnBtmDet','OffBtmDet'})
legend('location','eastoutside')
title('Bottom Detritus Flux (gC m^-^2 d^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_mean_feisty_forcing.png'])

%% Vert
figure(3)
subplot(2,2,3)
plot(log10(NvMZ(:,1)+eps),-1*z_l,'color',cm10(4,:)); hold on;
plot(log10(NvLZ(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on;
plot(log10(FvMZ(:,1)+eps),-1*z_l,'--','color',cm10(4,:)); hold on;
plot(log10(FvLZ(:,1)+eps),-1*z_l,'--','color',cm10(5,:)); hold on;
title('log_1_0')
%legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
subplot(2,2,4)
plot(log10(NvMH(:,1)+eps),-1*z_l,'color',cm10(6,:)); hold on; 
plot(log10(NvLH(:,1)+eps),-1*z_l,'color',cm10(7,:)); hold on;
plot(log10(FvMH(:,1)+eps),-1*z_l,'--','color',cm10(6,:)); hold on; 
plot(log10(FvLH(:,1)+eps),-1*z_l,'--','color',cm10(7,:)); hold on;
title('log_1_0')
%legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
%legend('location','east')
ylabel('Depth (m)')

subplot(2,2,1)
plot((NvMZ(1:8,1)),-1*z_l(1:8),'color',cm10(4,:)); hold on;
plot((NvLZ(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on;
plot((FvMZ(1:8,1)),-1*z_l(1:8),'--','color',cm10(4,:)); hold on;
plot((FvLZ(1:8,1)),-1*z_l(1:8),'--','color',cm10(5,:)); hold on;
legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
subplot(2,2,2)
plot((NvMH(1:8,1)),-1*z_l(1:8),'color',cm10(6,:)); hold on; 
plot((NvLH(1:8,1)),-1*z_l(1:8),'color',cm10(7,:)); hold on;
plot((FvMH(1:8,1)),-1*z_l(1:8),'--','color',cm10(6,:)); hold on; 
plot((FvLH(1:8,1)),-1*z_l(1:8),'--','color',cm10(7,:)); hold on;
legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
legend('location','southeast')
title('Mean Higher Predation Rate (gC m^-^3 d^-^1)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_mean_feisty_forcing.png'])

%% Maps
% Zoo bio 
figure(4)
subplot(2,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsMZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean Medium zoo (gC m^-^2)','HorizontalAlignment','center')

subplot(2,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsMZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')

subplot(2,2,3) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsLZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean Large zoo (gC m^-^2)','HorizontalAlignment','center')

subplot(2,2,4) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsLZ(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_zmeso.png'])


% Zoo HPloss
figure(5)
subplot(2,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsMH(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean MZ HPloss (gC m^-^2 d^-^1)','HorizontalAlignment','center')

subplot(2,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsMH(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')

subplot(2,2,3) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsLH(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(6,2.5,'log10 mean LZ HPloss (gC m^-^2 d^-^1)','HorizontalAlignment','center')

subplot(2,2,4) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsLH(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_hploss.png'])



% Det btm
figure(6)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(FsDE(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Offline')
text(4,2.5,'log10 mean btm detritus flux (gC m^-^2 d^-^1)','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(NsBDE(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Online')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_DetBtm.png'])


%% PLOTS OF DIFFS
% Time series
% NO Log10
figure(7)
subplot(3,1,1)
plot(tmos,(dtMZ(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on;
plot(tmos,(dtLZ(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
title('Integrated Biomass (gC m^-^2)')
legend({'MZ','LZ'})
legend('location','east')
subplot(3,1,2)
plot(tmos,(dtMH(tmos)* 365+eps),'color',cm10(6,:),'LineWidth',2); hold on; 
plot(tmos,(dtLH(tmos)* 365+eps),'color',cm10(7,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated Higher Predation Rate (gC m^-^2 y^-^1)')
legend({'MZ','LZ'})
legend('location','northeast')
subplot(3,1,3)
plot(tmos,(dtDE(tmos)* 365+eps),'color',cm10(1,:),'LineWidth',2); 
title('Bottom Detritus Flux (gC m^-^2 y^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_diff_feisty_forcing.png'])

%% Vert
figure(8)
% subplot(2,2,3)
% plot((dvMZ(:,1)+eps),-1*z_l,'color',cm10(4,:)); hold on;
% plot((dvLZ(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on;
% title('Online - Offline')
% subplot(2,2,4)
% plot((dvMH(:,1)+eps),-1*z_l,'color',cm10(6,:)); hold on; 
% plot((dvLH(:,1)+eps),-1*z_l,'color',cm10(7,:)); hold on;
% title('Online - Offline')
% ylabel('Depth (m)')

subplot(1,2,1)
plot((dvMZ(1:8,1)),-1*z_l(1:8),'color',cm10(4,:),'LineWidth',2); hold on;
plot((dvLZ(1:8,1)),-1*z_l(1:8),'color',cm10(5,:),'LineWidth',2); hold on;
legend({'MZ','LZ'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
subplot(1,2,2)
plot((dvMH(1:8,1)* 365),-1*z_l(1:8),'color',cm10(6,:),'LineWidth',2); hold on; 
plot((dvLH(1:8,1)* 365),-1*z_l(1:8),'color',cm10(7,:),'LineWidth',2); hold on;
legend({'MZ','LZ'})
legend('location','southwest')
title('Mean Higher Predation Rate (gC m^-^3 y^-^1)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_diff_feisty_forcing.png'])


%% Maps
% MZ & LZ
figure(9)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsMZ(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Medium zoo biomass difference (gC m^2)')
text(4,2.5,'Online - Offline','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLZ(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Large zoo biomass difference (gC m^2)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_zmeso.png'])

%% HPloss
figure(10)
subplot(1,2,1) %off
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsMH(:,:,1))* 365))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('MZ HPloss difference (gC m^2 y^-^1)')
text(4,2.5,'Online - Offline','HorizontalAlignment','center')

subplot(1,2,2) %on
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLH(:,:,1))* 365))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('LZ HPloss difference (gC m^2 y^-^1)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_HPloss.png'])

%% Det
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDE(:,:,1))* 365))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-25 25]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Bottom detrius flux difference (gC m^2 y^-^1) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_DetBtm.png'])

