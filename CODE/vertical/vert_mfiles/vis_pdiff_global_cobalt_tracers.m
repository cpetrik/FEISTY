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

%% Percent Diffs Online - Offline / Offline

dtDE = 100* (NtDE - FtDE) ./ FtDE;
dtDI = 100* (NtDI - FtDI) ./ FtDI;
dtSP = 100* (NtSP - FtSP) ./ FtSP;
dtLP = 100* (NtLP - FtLP) ./ FtLP;
dtSZ = 100* (NtSZ - FtSZ) ./ FtSZ;

dsDE = 100* (NsDE - FsDE) ./ FsDE;
dsDI = 100* (NsDI - FsDI) ./ FsDI;
dsSP = 100* (NsSP - FsSP) ./ FsSP;
dsLP = 100* (NsLP - FsLP) ./ FsLP;
dsSZ = 100* (NsSZ - FsSZ) ./ FsSZ;

dvDE = 100* (NvDE - FvDE) ./ FvDE;
dvDI = 100* (NvDI - FvDI) ./ FvDI;
dvSP = 100* (NvSP - FvSP) ./ FvSP;
dvLP = 100* (NvLP - FvLP) ./ FvLP;
dvSZ = 100* (NvSZ - FvSZ) ./ FvSZ;

tmos = 1:60;

%% PLOTS OF PERCENT DIFFS ---------------------------------------
% Time series
% NO Log10
figure(7)
subplot(3,1,1)
plot(tmos,(dtSP(tmos)+eps),'color',cm10(3,:),'LineWidth',2); hold on;
plot(tmos,(dtLP(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on; 
plot(tmos,(dtDI(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
title('Integrated Biomass (gC m^-^2) %Difference Online - Offline')
legend({'SP','LP','Diaz'})
legend('location','southwest')
subplot(3,1,2)
plot(tmos,(dtSZ(tmos)+eps),'color',cm10(8,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated SZ Biomass (gC m^-^2) %Difference Online - Offline')
% xlabel('Months')
subplot(3,1,3)
plot(tmos,(dtDE(tmos)+eps),'color',cm10(1,:),'LineWidth',2); 
title('Detritus Concentration (gC m^-^2) %Difference Online - Offline')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_Pdiff_cobalt_tracers.png'])

%% Vert
figure(8)
subplot(1,2,1)
plot((dvSP(1:9,1)),-1*z_l(1:9),'color',cm10(3,:),'LineWidth',2); hold on;
plot((dvLP(1:9,1)),-1*z_l(1:9),'color',cm10(4,:),'LineWidth',2); hold on; 
plot((dvDI(1:9,1)),-1*z_l(1:9),'color',cm10(5,:),'LineWidth',2); hold on;
legend({'SP','LP','Diaz'})
legend('location','southwest')
title('Mean Biomass (gC m^-^3) %Difference Online - Offline')
subplot(1,2,2)
plot((dvSZ(1:9,1)),-1*z_l(1:9),'color',cm10(8,:),'LineWidth',2); hold on;
title('Mean Biomass (gC m^-^3) %Difference Online - Offline')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_Pdiff_cobalt_tracers.png'])


%% Maps
% SZ
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsSZ(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Small zoo %difference (gC m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_Pdiff_SZ.png'])

%% Det
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDE(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Detrius %difference (gC m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_Pdiff_Det.png'])

%% SP
figure(12) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsSP(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Small phyto biomass %difference (gC m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_SP.png'])

%% LP
figure(13)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLP(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Large phyto biomass %difference (gC m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_LP.png'])

%% DI
figure(14)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDI(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Diazotroph biomass %difference (gC m^2)')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_Diaz.png'])



