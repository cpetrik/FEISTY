% Plot differences between FEISTY only (2D) and COBALT-FEISTY (3D)
% Feisty biomasses

clear
close all

%% ONLINE -----------------------------------------------------------
npath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

load([npath '19900101.ocean_feisty_tracers_z_means.mat'])

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

FZsf=NaN*ones(ni,nj);
FZsp=NaN*ones(ni,nj);
FZsd=NaN*ones(ni,nj);
FZmf=NaN*ones(ni,nj);
FZmp=NaN*ones(ni,nj);
FZmd=NaN*ones(ni,nj);
FZlp=NaN*ones(ni,nj);
FFZld=NaN*ones(ni,nj);
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
FAll = FZsp+FZsf+FZsd+FZmp+FZmf+FZmd+FZlp+FZld;
FsF = Zsf+FZmf;
FsP = FZsp+FZmp+FZlp;
FsD = FZsd+FZmd+FZld;
FsM = FZmp+FZmf+FZmd;
FsL = FZlp+FZld;
FFracPD = FsP ./ (FsP + FsD);
FFracPF = FsP ./ (FsP + FsF);
FFracLM = FsL ./ (FsL + FsM);


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

%% PLOTS TOGETHER
% Time series
tmos = 1:12;

% Log10
figure(1)
subplot(3,1,1)
plot(tmos,log10(NtMZ(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on;
plot(tmos,log10(NtLZ(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
plot(tmos,log10(FtMZ(tmos)+eps),'--','color',cm10(4,:),'LineWidth',2); hold on;
plot(tmos,log10(FtLZ(tmos)+eps),'--','color',cm10(5,:),'LineWidth',2); hold on;
legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
legend('location','eastoutside')
title('log_1_0 Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
plot(tmos,log10(NtMH(tmos)+eps),'color',cm10(6,:),'LineWidth',2); hold on; 
plot(tmos,log10(NtLH(tmos)+eps),'color',cm10(7,:),'LineWidth',2); hold on;
plot(tmos,log10(FtMH(tmos)+eps),'--','color',cm10(6,:),'LineWidth',2); hold on; 
plot(tmos,log10(FtLH(tmos)+eps),'--','color',cm10(7,:),'LineWidth',2); hold on;
legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
legend('location','eastoutside')
title('log_1_0 Integrated Higher Predation Rate (gC m^-^2 d^-^1)')
subplot(3,1,3)
plot(tmos,log10(NtDE(tmos)+eps),'color',cm10(1,:),'LineWidth',2); hold on;
plot(tmos,log10(FtDE(tmos)+eps),'--','color',cm10(1,:),'LineWidth',2); hold on;
legend({'OnBtmDet','OffBtmDet'})
legend('location','eastoutside')
title('log_1_0 Bottom Detritus Flux (gC m^-^2 d^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_log10mean_feisty_forcing.png'])

% No log
figure(2)
subplot(3,1,1)
plot(tmos,(NtMZ(tmos)+eps),'color',cm10(4,:),'LineWidth',2); hold on;
plot(tmos,(NtLZ(tmos)+eps),'color',cm10(5,:),'LineWidth',2); hold on;
plot(tmos,(FtMZ(tmos)+eps),'--','color',cm10(4,:),'LineWidth',2); hold on;
plot(tmos,(FtLZ(tmos)+eps),'--','color',cm10(5,:),'LineWidth',2); hold on;
legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
legend('location','eastoutside')
title('Integrated Biomass (gC m^-^2)')
subplot(3,1,2)
plot(tmos,(NtMH(tmos)+eps),'color',cm10(6,:),'LineWidth',2); hold on; 
plot(tmos,(NtLH(tmos)+eps),'color',cm10(7,:),'LineWidth',2); hold on;
plot(tmos,(FtMH(tmos)+eps),'--','color',cm10(6,:),'LineWidth',2); hold on; 
plot(tmos,(FtLH(tmos)+eps),'--','color',cm10(7,:),'LineWidth',2); hold on;
legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
legend('location','eastoutside')
title('Integrated Higher Predation Rate (gC m^-^2 d^-^1)')
subplot(3,1,3)
plot(tmos,(NtDE(tmos)+eps),'color',cm10(1,:),'LineWidth',2); hold on;
plot(tmos,(FtDE(tmos)+eps),'--','color',cm10(1,:),'LineWidth',2); hold on;
legend({'OnBtmDet','OffBtmDet'})
legend('location','eastoutside')
title('Bottom Detritus Flux (gC m^-^2 d^-^1)')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_mean_feisty_forcing.png'])

%% Vert
figure(3)
subplot(2,2,3)
plot(log10(NvMZ(:,1)+eps),-1*z_l,'color',cm10(4,:),'LineWidth',2); hold on;
plot(log10(NvLZ(:,1)+eps),-1*z_l,'color',cm10(5,:),'LineWidth',2); hold on;
plot(log10(FvMZ(:,1)+eps),-1*z_l,'--','color',cm10(4,:),'LineWidth',2); hold on;
plot(log10(FvLZ(:,1)+eps),-1*z_l,'--','color',cm10(5,:),'LineWidth',2); hold on;
title('log_1_0')
%legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
subplot(2,2,4)
plot(log10(NvMH(:,1)+eps),-1*z_l,'color',cm10(6,:),'LineWidth',2); hold on; 
plot(log10(NvLH(:,1)+eps),-1*z_l,'color',cm10(7,:),'LineWidth',2); hold on;
plot(log10(FvMH(:,1)+eps),-1*z_l,'--','color',cm10(6,:),'LineWidth',2); hold on; 
plot(log10(FvLH(:,1)+eps),-1*z_l,'--','color',cm10(7,:),'LineWidth',2); hold on;
title('log_1_0')
%legend({'OnMZloss','OnLZloss','OffMZloss','OffLZloss'})
%legend('location','east')
ylabel('Depth (m)')

subplot(2,2,1)
plot((NvMZ(1:8,1)),-1*z_l(1:8),'color',cm10(4,:),'LineWidth',2); hold on;
plot((NvLZ(1:8,1)),-1*z_l(1:8),'color',cm10(5,:),'LineWidth',2); hold on;
plot((FvMZ(1:8,1)),-1*z_l(1:8),'--','color',cm10(4,:),'LineWidth',2); hold on;
plot((FvLZ(1:8,1)),-1*z_l(1:8),'--','color',cm10(5,:),'LineWidth',2); hold on;
legend({'OnMZbio','OnLZbio','OffMZbio','OffLZbio'})
legend('location','southeast')
title('Mean Biomass (gC m^-^3)')
subplot(2,2,2)
plot((NvMH(1:8,1)),-1*z_l(1:8),'color',cm10(6,:),'LineWidth',2); hold on; 
plot((NvLH(1:8,1)),-1*z_l(1:8),'color',cm10(7,:),'LineWidth',2); hold on;
plot((FvMH(1:8,1)),-1*z_l(1:8),'--','color',cm10(6,:),'LineWidth',2); hold on; 
plot((FvLH(1:8,1)),-1*z_l(1:8),'--','color',cm10(7,:),'LineWidth',2); hold on;
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
surfm(geolat,geolon,log10(squeeze(NsDE(:,:,1))))
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

