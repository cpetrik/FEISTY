% Plot differences between COBALT only and COBALT-FEISTY 
% COBALT nuts & phyto vars

clear
close all

%% molN/kg to gC/m3
NtoC= 1035 * (106/16) * 12.01;

%% ONLINE -----------------------------------------------------------
npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

load([npath 'ocean_cobalt_tracers_month_z.199001-199412_means.mat'])

%% Put everything in gC m^-^2 
NtDE = NtoC * tDE(1:60);
NtDI = NtoC * tDI(1:60);
NtSP = NtoC * tSP(1:60);
NtLP = NtoC * tLP(1:60);
NtSZ = NtoC * tSZ(1:60);

NvDE = NtoC * vDE(:,1);
NvDI = NtoC * vDI(:,1);
NvSP = NtoC * vSP(:,1);
NvLP = NtoC * vLP(:,1);
NvSZ = NtoC * vSZ(:,1);

NsDE = NtoC * sDE(:,:,1);
NsDI = NtoC * sDI(:,:,1);
NsSP = NtoC * sSP(:,:,1);
NsLP = NtoC * sLP(:,:,1);
NsSZ = NtoC * sSZ(:,:,1);

clear tDE tDI vDE vDI sDE sDI
clear tSP tLP tSZ vSP vLP vSZ sSP sLP sSZ

%% OFFLINE -----------------------------------------------------------
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_tracers_month_z.199001-199412_means.mat'])

%% Put everything in gC m^-^2 
% FtDE = NtoC * tDE(1:60);
% FtDI = NtoC * tDI(1:60);
FtSP = NtoC * tSP(1:60);
FtLP = NtoC * tLP(1:60);
FtSZ = NtoC * tSZ(1:60);

% FvDE = NtoC * vDE(:,1);
% FvDI = NtoC * vDI(:,1);
FvSP = NtoC * vSP(:,1);
FvLP = NtoC * vLP(:,1);
FvSZ = NtoC * vSZ(:,1);

% FsDE = NtoC * sDE(:,:,1);
% FsDI = NtoC * sDI(:,:,1);
FsSP = NtoC * sSP(:,:,1);
FsLP = NtoC * sLP(:,:,1);
FsSZ = NtoC * sSZ(:,:,1);

clear tDE tDI vDE vDI sDE sDI
clear tSP tLP tSZ vSP vLP vSZ sSP sLP sSZ

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

% dtDE = NtDE - FtDE;
% dtDI = NtDI - FtDI;
dtSP = NtSP - FtSP;
dtLP = NtLP - FtLP;
dtSZ = NtSZ - FtSZ;

% dsDE = NsDE - FsDE;
% dsDI = NsDI - FsDI;
dsSP = NsSP - FsSP;
dsLP = NsLP - FsLP;
dsSZ = NsSZ - FsSZ;

% dvDE = NvDE - FvDE;
% dvDI = NvDI - FvDI;
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
plot(tmos,log10(FtSP(tmos)+eps),'--','color',cm10(3,:)); hold on;
plot(tmos,log10(FtLP(tmos)+eps),'--','color',cm10(4,:)); hold on; 
legend({'OnSP','OnLP','OffSP','OffLP'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
%plot(tmos,log10(NtDI(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,log10(NtSZ(tmos)+eps),'color',cm10(8,:)); hold on;
%plot(tmos,log10(FtDI(tmos)+eps),'--','color',cm10(5,:)); hold on;
plot(tmos,log10(FtSZ(tmos)+eps),'--','color',cm10(8,:)); hold on;
%legend({'OnDI','OnSZ','OffDI','OffSZ'})
legend({'OnSZ','OffSZ'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^2)')
xlabel('Months')
% subplot(3,1,3)
% plot(tmos,log10(NtDE(tmos)+eps),'color',cm10(1,:)); hold on;
% plot(tmos,log10(FtDE(tmos)+eps),'--','color',cm10(1,:)); hold on;
% legend({'OnDet','OffDet'})
% legend('location','eastoutside')
% title('log_1_0 Integrated Detritus (gC m^-^2)')
% xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_log10mean_cobalt_tracers.png'])

%% No log
figure(2)
subplot(3,1,1)
plot(tmos,(NtSP(tmos)+eps),'color',cm10(3,:)); hold on;
plot(tmos,(NtLP(tmos)+eps),'color',cm10(4,:)); hold on; 
plot(tmos,(FtSP(tmos)+eps),'--','color',cm10(3,:)); hold on;
plot(tmos,(FtLP(tmos)+eps),'--','color',cm10(4,:)); hold on; 
legend({'OnSP','OnLP','OffSP','OffLP'})
legend('location','eastoutside')
title('Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
%plot(tmos,(NtDI(tmos)+eps),'color',cm10(5,:)); hold on;
plot(tmos,(NtSZ(tmos)+eps),'color',cm10(8,:)); hold on;
%plot(tmos,(FtDI(tmos)+eps),'--','color',cm10(5,:)); hold on;
plot(tmos,(FtSZ(tmos)+eps),'--','color',cm10(8,:)); hold on;
%legend({'OnDI','OnSZ','OffDI','OffSZ'})
legend({'OnSZ','OffSZ'})
legend('location','eastoutside')
title('Integrated Biomass (gC m^-^2)')
xlabel('Months')
% subplot(3,1,3)
% plot(tmos,(NtDE(tmos)+eps),'color',cm10(1,:)); hold on;
% plot(tmos,(FtDE(tmos)+eps),'--','color',cm10(1,:)); hold on;
% legend({'OnDet','OffDet'})
% legend('location','eastoutside')
% title('Integrated Detritus (gC m^-^2)')
% xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_mean_cobalt_tracers.png'])

%% Vert
figure(3)
subplot(2,2,3)
plot(log10(NvSP(:,1)+eps),-1*z_l,'color',cm10(3,:)); hold on;
plot(log10(NvLP(:,1)+eps),-1*z_l,'color',cm10(4,:)); hold on; 
plot(log10(FvSP(:,1)+eps),-1*z_l,'--','color',cm10(3,:)); hold on;
plot(log10(FvLP(:,1)+eps),-1*z_l,'--','color',cm10(4,:)); hold on; 
title('log_1_0')
subplot(2,2,4)
% plot(log10(NvDI(:,1)+eps),-1*z_l,'color',cm10(5,:)); hold on;
plot(log10(NvSZ(:,1)+eps),-1*z_l,'color',cm10(8,:)); hold on;
% plot(log10(FvDI(:,1)+eps),-1*z_l,'--','color',cm10(5,:)); hold on;
plot(log10(FvSZ(:,1)+eps),-1*z_l,'--','color',cm10(8,:)); hold on;
title('log_1_0')
ylabel('Depth (m)')

subplot(2,2,1)
plot((NvSP(1:8,1)),-1*z_l(1:8),'color',cm10(3,:)); hold on;
plot((NvLP(1:8,1)),-1*z_l(1:8),'color',cm10(4,:)); hold on; 
plot((FvSP(1:8,1)),-1*z_l(1:8),'--','color',cm10(3,:)); hold on;
plot((FvLP(1:8,1)),-1*z_l(1:8),'--','color',cm10(4,:)); hold on; 
legend({'OnSP','OnLP','OffSP','OffLP'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
subplot(2,2,2)
% plot((NvDI(1:8,1)),-1*z_l(1:8),'color',cm10(5,:)); hold on;
plot((NvSZ(1:8,1)),-1*z_l(1:8),'color',cm10(8,:)); hold on;
% plot((FvDI(1:8,1)),-1*z_l(1:8),'--','color',cm10(5,:)); hold on;
plot((FvSZ(1:8,1)),-1*z_l(1:8),'--','color',cm10(8,:)); hold on;
% legend({'OnDIs','OnSZ','OffDI','OffSZ'})
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
% figure(6)
% subplot(1,2,1) %off
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,log10(squeeze(FsDI(:,:,1))))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-1 1]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Offline')
% text(4,2.5,'log10 mean Diaz (gC m^-^2)','HorizontalAlignment','center')
% 
% subplot(1,2,2) %on
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,log10(squeeze(NsDI(:,:,1))))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-1 1]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Online')
% stamp('')
% print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Diaz.png'])


%% Det 
% figure(7)
% subplot(1,2,1) %off
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,log10(squeeze(FsDE(:,:,1))))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-3 -1]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Offline')
% text(4,2.5,'log10 mean btm detritus flux (gC m^-^2 d^-^1)','HorizontalAlignment','center')
% 
% subplot(1,2,2) %on
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,log10(squeeze(NsDE(:,:,1))))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-3 -1]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Online')
% stamp('')
% print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_Det.png'])


%% PLOTS OF DIFFS ---------------------------------------
% Time series
% NO Log10
figure(7)
subplot(3,1,1)
plot(tmos,(dtSP(tmos)+eps),'color',cm10(3,:),'LineWidth',2); hold on;
plot(tmos,(dtLP(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on; 
title('Integrated Biomass (gC m^-^2) Difference Online - Offline')
legend({'SP','LP'})
legend('location','northeast')
subplot(3,1,2)
%plot(tmos,(dtDI(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
plot(tmos,(dtSZ(tmos)+eps),'color',cm10(8,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated Biomass (gC m^-^2) Difference Online - Offline')
%legend({'DI','SZ'})
%legend('location','northeast')
xlabel('Months')
% subplot(3,1,3)
% plot(tmos,(dtDE(tmos)+eps),'color',cm10(1,:),'LineWidth',2); 
% title('Detritus Concentration (gC m^-^2) Difference Online - Offline')
% xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_diff_cobalt_tracers.png'])

%% Vert
figure(8)
subplot(1,2,1)
plot((dvSP(1:9,1)),-1*z_l(1:9),'color',cm10(3,:),'LineWidth',2); hold on;
plot((dvLP(1:9,1)),-1*z_l(1:9),'color',cm10(4,:),'LineWidth',2); hold on; 
legend({'SP','LP'})
legend('location','southwest')
title('Mean Biomass (gC m^-^3) Difference Online - Offline')
subplot(1,2,2)
% plot((dvDI(1:9,1)),-1*z_l(1:9),'color',cm10(5,:),'LineWidth',2); hold on;
plot((dvSZ(1:9,1)),-1*z_l(1:9),'color',cm10(8,:),'LineWidth',2); hold on;
% legend({'DI','SZ'})
% legend('location','southwest')
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

% subplot(2,2,3) %DI
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,(squeeze(dsDI(:,:,1))))
% cmocean('balance')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-3 3]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Diazotroph biomass difference (gC m^2)')
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
% figure(11)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,(squeeze(dsDE(:,:,1))))
% cmocean('balance')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-5 5]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Detrius difference (gC m^2) Online - Offline')
% stamp('')
% print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_Det.png'])

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
%figure(14)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat,geolon,(squeeze(dsDI(:,:,1))))
% cmocean('balance')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-3 3]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Diazotroph biomass difference (gC m^2)')
%print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_Diaz.png'])



