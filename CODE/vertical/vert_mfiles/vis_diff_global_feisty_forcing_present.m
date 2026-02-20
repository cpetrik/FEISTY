% Plot differences between COBALT only and COBALT-FEISTY 
% Feisty forcing vars

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

%% ONLINE -----------------------------------------------------------
npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

load([npath 'ocean_cobalt_feisty_forcing_z.199001-199912_means.mat'])

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
NtDE = N2CmY * tDE;
NtMZ = N2Ckg * tMZ;
NtLZ = N2Ckg * tLZ;
NtMH = N2CkgY * tMH;
NtLH = N2CkgY * tLH;

NsDE = N2CmY * sDE;
NsMZ = N2Ckg * sMZ;
NsLZ = N2Ckg * sLZ;
NsMH = N2CkgY * sMH;
NsLH = N2CkgY * sLH;

NvMZ = N2Ckg * vMZ;
NvLZ = N2Ckg * vLZ;
NvMH = N2CkgY * vMH;
NvLH = N2CkgY * vLH;

clear tDE tMZ tLZ tMH tLH
clear sDE sMZ sLZ sMH sLH
clear vMZ vLZ vMH vLH

%% OFFLINE -----------------------------------------------------------
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_feisty_forcing_z.199001-199912_means.mat'])

%% Put everything in gC m^-^2 or gC m^-^2 d^-^1
FtDE = N2CmY * tDE;
FtMZ = N2Ckg * tMZ;
FtLZ = N2Ckg * tLZ;
FtMH = N2CkgY * tMH;
FtLH = N2CkgY * tLH;

FsDE = N2CmY * sDE;
FsMZ = N2Ckg * sMZ;
FsLZ = N2Ckg * sLZ;
FsMH = N2CkgY * sMH;
FsLH = N2CkgY * sLH;

FvMZ = N2Ckg * vMZ;
FvLZ = N2Ckg * vLZ;
FvMH = N2CkgY * vMH;
FvLH = N2CkgY * vLH;

clear tDE tMZ tLZ tMH tLH
clear sDE sMZ sLZ sMH sLH
clear vMZ vLZ vMH vLH

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

dtDE = NtDE - FtDE;
dtMZ = NtMZ - FtMZ;
dtLZ = NtLZ - FtLZ;
dtMH = NtMH - FtMH;
dtLH = NtLH - FtLH;

dsDE = NsDE - FsDE;
dsMZ = NsMZ - FsMZ;
dsLZ = NsLZ - FsLZ;
dsMH = NsMH - FsMH;
dsLH = NsLH - FsLH;

dvMZ = NvMZ - FvMZ;
dvLZ = NvLZ - FvLZ;
dvMH = NvMH - FvMH;
dvLH = NvLH - FvLH;


%% PLOTS OF DIFFS
% Time series
% NO Log10
figure(1)
subplot(3,1,1)
plot(tmos,(dtMZ(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on;
plot(tmos,(dtLZ(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
title('Integrated Biomass (gC m^-^2)')
legend({'MZ','LZ'})
legend('location','east')
subplot(3,1,2)
plot(tmos,(dtMH(tmos)+eps),'color',cm10(6,:),'LineWidth',2); hold on; 
plot(tmos,(dtLH(tmos)+eps),'color',cm10(7,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated Higher Predation Rate (gC m^-^2 y^-^1)')
legend({'MZ','LZ'})
legend('location','northeast')
subplot(3,1,3)
plot(tmos,(dtDE(tmos)+eps),'color',cm10(1,:),'LineWidth',2); 
title('Bottom Detritus Flux (gC m^-^2 y^-^1)')
xlabel('Months')
%stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_diff_feisty_forcing.png'])

%% Vert
figure(2)
subplot(1,2,1)
plot((dvMZ(1:8,1)),-1*z_l(1:8),'color',cm10(4,:),'LineWidth',2); hold on;
plot((dvLZ(1:8,1)),-1*z_l(1:8),'color',cm10(5,:),'LineWidth',2); hold on;
legend({'MZ','LZ'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
subplot(1,2,2)
plot((dvMH(1:8,1)),-1*z_l(1:8),'color',cm10(6,:),'LineWidth',2); hold on; 
plot((dvLH(1:8,1)),-1*z_l(1:8),'color',cm10(7,:),'LineWidth',2); hold on;
legend({'MZ','LZ'})
legend('location','southwest')
title('Mean Higher Predation Rate (gC m^-^3 y^-^1)')
ylabel('Depth (m)')
%stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_diff_feisty_forcing.png'])

%% Maps
% MZ & LZ
figure(3) %MZ
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsMZ(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Medium zoo biomass difference (gC m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_MZ.png'])

figure(4) %LZ
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLZ(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-3 3]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Large zoo biomass difference (gC m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_diff_LZ.png'])

%% HPloss
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsMH(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('MZ HPloss difference (gC m^2 y^-^1) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_mzHPloss.png'])

figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsLH(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('LZ HPloss difference (gC m^2 y^-^1) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_lzHPloss.png'])

%% Det
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsDE(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-25 25]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Bottom detrius flux difference (gC m^2 y^-^1) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_diff_DetBtm.png'])

