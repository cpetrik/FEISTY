% Look at COBALT-FEISTY fluxes (fish) from online sim
% consump zoop terms

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';
mod = exper;
Yr = '1990'; %1990, 1994

%%
load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet',...
    'z_l_units','z_l_long_name','z_l','geolon','geolat')

dz = diff(z_l);
%dz_mat = repmat(dz,1,12);

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
LON = double(geolon);
LAT = double(geolat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% 
load([fpath Yr '0101.ocean_feisty_consump_z_means.mat'])

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
%y = 1:12; %120;

%% Time series
mos = 1:12;

% MZ
figure(1)
subplot(3,1,1)
plot(mos,(ts_means.SF_cons_Mz),'r'); hold on;
plot(mos,(ts_means.SP_cons_Mz),'--b'); hold on;
plot(mos,(ts_means.SD_cons_Mz),':g'); hold on;
legend({'SF','SP','SD'})
title('Mean MZ consump (d^-^1)')

subplot(3,1,2)
plot(mos,(ts_means.MF_cons_Mz),'r'); hold on;
plot(mos,(ts_means.MP_cons_Mz),'--b'); hold on;
plot(mos,(ts_means.MD_cons_Mz),':g'); hold on;
legend({'MF','MP','MD'})

subplot(3,1,3)
plot(mos,(ts_means.LP_cons_Mz),'b'); hold on;
plot(mos,(ts_means.LD_cons_Mz),'--g'); hold on;
legend({'LP','LD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_' Yr '_cons_Mz_ts_mean_fntypes.png'])

%% LZ consump
figure(2)
subplot(3,1,1)
plot(mos,(ts_means.SF_cons_Lz),'r'); hold on;
plot(mos,(ts_means.SP_cons_Lz),'--b'); hold on;
plot(mos,(ts_means.SD_cons_Lz),':g'); hold on;
legend({'SF','SP','SD'})
title('Mean LZ consump (d^-^1)')

subplot(3,1,2)
plot(mos,(ts_means.MF_cons_Lz),'r'); hold on;
plot(mos,(ts_means.MP_cons_Lz),'--b'); hold on;
plot(mos,(ts_means.MD_cons_Lz),':g'); hold on;
legend({'MF','MP','MD'})

subplot(3,1,3)
plot(mos,(ts_means.LP_cons_Lz),'b'); hold on;
plot(mos,(ts_means.LD_cons_Lz),'--g'); hold on;
legend({'LP','LD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_' Yr '_cons_Lz_ts_mean_fntypes.png'])

%% Bent consump
figure(3)
subplot(3,1,1)
plot(mos,(ts_means.SF_cons_BE),'r'); hold on;
plot(mos,(ts_means.SP_cons_BE),'--b'); hold on;
plot(mos,(ts_means.SD_cons_BE),':g'); hold on;
legend({'SF','SP','SD'})
title('Mean LZ consump (d^-^1)')

subplot(3,1,2)
plot(mos,(ts_means.MF_cons_BE),'r'); hold on;
plot(mos,(ts_means.MP_cons_BE),'--b'); hold on;
plot(mos,(ts_means.MD_cons_BE),':g'); hold on;
legend({'MF','MP','MD'})

subplot(3,1,3)
plot(mos,(ts_means.LP_cons_BE),'b'); hold on;
plot(mos,(ts_means.LD_cons_BE),'--g'); hold on;
legend({'LP','LD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_' Yr '_cons_BE_ts_mean_fntypes.png'])


%% Vert distrib
figure(6) % MZ
subplot(2,3,1)
plot((vert_means.SF_cons_Mz(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot((vert_means.MF_cons_Mz(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
title('Forage')
ylabel('Depth (m)')

subplot(2,3,2)
plot((vert_means.SP_cons_Mz(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot((vert_means.MP_cons_Mz(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot((vert_means.LP_cons_Mz(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
title('Lg pel')
xlabel('MZ consump (d^-^1)')

subplot(2,3,3)
plot((vert_means.SD_cons_Mz(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
title('Demersal')

subplot(2,3,4)
plot(log10(vert_means.SF_cons_Mz),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vert_means.MF_cons_Mz),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
legend('SF','MF')
%legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
xlabel('log_1_0 MZ consump (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,5)
plot(log10(vert_means.SP_cons_Mz),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vert_means.MP_cons_Mz),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vert_means.LP_cons_Mz),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
legend('SP','MP','LP')
legend('location','east')

subplot(2,3,6)
plot(log10(vert_means.SD_cons_Mz),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
plot(log10(vert_means.MD_cons_Mz),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
plot(log10(vert_means.LD_cons_Mz),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
legend('SD','MD','LD')
legend('location','east')
stamp('')
print('-dpng',[ppath mod '_' Yr '_cons_Mz_vert_mean_fntypes.png'])

%% vert LZ
figure(7) 
subplot(2,3,1)
plot((vert_means.SF_cons_Lz(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot((vert_means.MF_cons_Lz(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
title('Forage')
xlabel('LZ consump (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,2)
plot((vert_means.SP_cons_Lz(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot((vert_means.MP_cons_Lz(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot((vert_means.LP_cons_Lz(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
title('Lg pel')

subplot(2,3,3)
plot((vert_means.SD_cons_Lz(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
title('Demersal')

subplot(2,3,4)
plot(log10(vert_means.SF_cons_Lz),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vert_means.MF_cons_Lz),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
legend('SF','MF')
%legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
xlabel('log_1_0 LZ consump (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,5)
plot(log10(vert_means.SP_cons_Lz),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vert_means.MP_cons_Lz),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vert_means.LP_cons_Lz),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
legend('SP','MP','LP')
legend('location','east')

subplot(2,3,6)
plot(log10(vert_means.SD_cons_Lz),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
plot(log10(vert_means.MD_cons_Lz),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
plot(log10(vert_means.LD_cons_Lz),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
legend('SD','MD','LD')
legend('location','east')
stamp('')
print('-dpng',[ppath mod '_' Yr '_cons_Lz_vert_mean_fntypes.png'])

%% vert BE
figure(8) 
subplot(2,3,1)
plot((vert_means.SF_cons_BE(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot((vert_means.MF_cons_BE(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
title('Forage')

subplot(2,3,2)
plot((vert_means.SP_cons_BE(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot((vert_means.MP_cons_BE(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot((vert_means.LP_cons_BE(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
title('Lg pel')

subplot(2,3,3)
plot((vert_means.SD_cons_BE(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
title('Demersal')

subplot(2,3,4)
plot(log10(vert_means.SF_cons_BE),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vert_means.MF_cons_BE),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
legend('SF','MF')
%legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
xlabel('log_1_0 ')
ylabel('Depth (m)')

subplot(2,3,5)
plot(log10(vert_means.SP_cons_BE),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vert_means.MP_cons_BE),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vert_means.LP_cons_BE),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
legend('SP','MP','LP')
legend('location','east')

subplot(2,3,6)
plot(log10(vert_means.SD_cons_BE),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
plot(log10(vert_means.MD_cons_BE),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
plot(log10(vert_means.LD_cons_BE),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
legend('SD','MD','LD')
legend('location','east')
stamp('')
%print('-dpng',[ppath mod '_' Yr '_cons_BE_vert_mean_fntypes.png'])

%% Maps

% All small on MZ
figure(9)
% SF
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.SF_cons_Mz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 MZ con SF (g m^-^2 d^-^1)')

% SD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.SD_cons_Mz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log10 MZ con SD (g m^-^2 d^-^1)')

% SP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.SP_cons_Mz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log10 MZ con SP (g m^-^2 d^-^1)')

% All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(LAT,LON,log10(sAll))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-2 2]);
% set(gcf,'renderer','painters')
% title('log10 mean All fishes (g m^-^2)')
stamp('')
print('-dpng',[ppath exper '_' Yr '_global_small_cons_Mz_subplot.png'])

%% All med MZ
figure(10)
% MF
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MF_cons_Mz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 MZ con MF (g m^-^2 d^-^1)')

% MD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MD_cons_Mz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log10 MZ con MD (g m^-^2 d^-^1)')

% MP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MP_cons_Mz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log10 MZ con MP (g m^-^2 d^-^1)')

% All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(LAT,LON,log10(sAll))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-2 2]);
% set(gcf,'renderer','painters')
% title('log10 mean All fishes (g m^-^2)')
stamp('')
print('-dpng',[ppath exper '_' Yr '_global_med_cons_Mz_subplot.png'])

%% All med LZ
figure(11)
% MF
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MF_cons_Lz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 LZ con MF (g m^-^2 d^-^1)')

% MD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MD_cons_Lz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log10 LZ con MD (g m^-^2 d^-^1)')

% MP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MP_cons_Lz))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log10 LZ con MP (g m^-^2 d^-^1)')

% All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(LAT,LON,log10(sAll))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-2 2]);
% set(gcf,'renderer','painters')
% title('log10 mean All fishes (g m^-^2)')
stamp('')
print('-dpng',[ppath exper '_' Yr '_global_med_cons_Lz_subplot.png'])


