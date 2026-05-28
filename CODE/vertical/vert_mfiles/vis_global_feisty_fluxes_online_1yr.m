% Look at COBALT-FEISTY fluxes (fish) from online sim

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
load([fpath '19940101.ocean_feisty_fluxes_z_means.mat'])

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
%y = 1:12; %120;

%% Time series
mos = 1:12;

% Met
figure(1)
subplot(3,1,1)
plot(mos,(ts_means.SF_met),'r'); hold on;
plot(mos,(ts_means.SP_met),'--b'); hold on;
plot(mos,(ts_means.SD_met),':g'); hold on;
legend({'SF','SP','SD'})
title('Mean Metab (d^-^1)')

subplot(3,1,2)
plot(mos,(ts_means.MF_met),'r'); hold on;
plot(mos,(ts_means.MP_met),'--b'); hold on;
plot(mos,(ts_means.MD_met),':g'); hold on;
legend({'MF','MP','MD'})

subplot(3,1,3)
plot(mos,(ts_means.LP_met),'b'); hold on;
plot(mos,(ts_means.LD_met),'--g'); hold on;
legend({'LP','LD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_met_ts_mean_fntypes.png'])

%% Rec
figure(2)
subplot(2,1,1)
plot(mos,(ts_means.SF_Fout),'r'); hold on;
plot(mos,(ts_means.SP_Fout),'--b'); hold on;
plot(mos,(ts_means.SD_Fout),':g'); hold on;
legend({'SF','SP','SD'})
title('Mean recruitment (d^-^1)')

subplot(2,1,2)
plot(mos,(ts_means.MP_Fout),'b'); hold on;
plot(mos,(ts_means.MD_Fout),'--g'); hold on;
legend({'MP','MD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_rec_ts_mean_fntypes.png'])

%% Prod
figure(3)
subplot(3,1,1)
plot(mos,(ts_means.SF_prod),'r'); hold on;
plot(mos,(ts_means.SP_prod),'--b'); hold on;
plot(mos,(ts_means.SD_prod),':g'); hold on;
legend({'SF','SP','SD'})
title('Mean production (gWW/m3/d)')

subplot(3,1,2)
plot(mos,(ts_means.MF_prod),'r'); hold on;
plot(mos,(ts_means.MP_prod),'--b'); hold on;
plot(mos,(ts_means.MD_prod),':g'); hold on;
legend({'MF','MP','MD'})

subplot(3,1,3)
plot(mos,(ts_means.LP_prod),'b'); hold on;
plot(mos,(ts_means.LD_prod),'--g'); hold on;
legend({'LP','LD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_prod_ts_mean_fntypes.png'])

%% Repro
figure(4)
subplot(2,1,1)
plot(mos,(ts_means.MF_rho),'r'); hold on;
title('Repro')
ylabel('MF (d^-^1)')

subplot(2,1,2)
plot(mos,(ts_means.LP_rho),'b'); hold on;
plot(mos,(ts_means.LD_rho),'--g'); hold on;
legend({'LP','LD'})
legend('location','southeast')
xlabel('Time (d)')
print('-dpng',[ppath exper '_repro_ts_mean_fntypes.png'])

%% mort
% figure(5)
% subplot(2,1,1)
% plot(mos,(ts_means.SF_mu),'r'); hold on;
% plot(mos,(ts_means.SP_mu),'--b'); hold on;
% plot(mos,(ts_means.SD_mu),':g'); hold on;
% legend({'SF','SP','SD'})
% title('Mean Predation mortality rate (d^-^1)')
% 
% subplot(2,1,2)
% plot(mos,(ts_means.MF_mu),'r'); hold on;
% plot(mos,(ts_means.MP_mu),'--b'); hold on;
% plot(mos,(ts_means.MD_mu),':g'); hold on;
% legend({'MF','MP','MD'})
% legend('location','southeast')
% xlabel('Time (d)')
% print('-dpng',[ppath exper '_pmort_ts_mean_fntypes.png'])


%% Vert distrib
figure(6) % Met
subplot(2,3,1)
plot((vert_means.SF_met(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot((vert_means.MF_met(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
title('Forage')
xlabel('Met (d^-^1)')

subplot(2,3,2)
plot((vert_means.SP_met(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot((vert_means.MP_met(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot((vert_means.LP_met(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
title('Lg pel')
xlabel('Met (d^-^1)')

subplot(2,3,3)
plot((vert_means.SD_met(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
title('Demersal')
xlabel('Met (d^-^1)')

subplot(2,3,4)
plot(log10(vert_means.SF_met),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vert_means.MF_met),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
legend('SF','MF')
%legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
xlabel('log_1_0 Met (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,5)
plot(log10(vert_means.SP_met),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vert_means.MP_met),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vert_means.LP_met),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
legend('SP','MP','LP')
legend('location','east')
xlabel('log_1_0 Met (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,6)
plot(log10(vert_means.SD_met),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
plot(log10(vert_means.MD_met),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
plot(log10(vert_means.LD_met),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
legend('SD','MD','LD')
legend('location','east')
xlabel('log_1_0 Met (d^-^1)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath mod '_met_vert_mean_fntypes.png'])

%% vert Rec
figure(7) 
subplot(2,3,1)
plot((vert_means.SF_Fout(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot((vert_means.MF_Fout(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
title('Forage')
xlabel('Rec (d^-^1)')

subplot(2,3,2)
plot((vert_means.SP_Fout(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot((vert_means.MP_Fout(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot((vert_means.LP_Fout(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
title('Lg pel')
xlabel('Rec (d^-^1)')

subplot(2,3,3)
plot((vert_means.SD_Fout(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
title('Demersal')
xlabel('Rec (d^-^1)')

subplot(2,3,4)
plot(log10(vert_means.SF_Fout),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vert_means.MF_Fout),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
legend('SF','MF')
%legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
xlabel('log_1_0 Rec (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,5)
plot(log10(vert_means.SP_Fout),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vert_means.MP_Fout),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vert_means.LP_Fout),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
legend('SP','MP','LP')
legend('location','east')
xlabel('log_1_0 Rec (d^-^1)')
ylabel('Depth (m)')

subplot(2,3,6)
plot(log10(vert_means.SD_Fout),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
plot(log10(vert_means.MD_Fout),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
plot(log10(vert_means.LD_Fout),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
legend('SD','MD','LD')
legend('location','east')
xlabel('log_1_0 Rec (d^-^1)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath mod '_recruit_vert_mean_fntypes.png'])

%% vert  Prod
figure(8) 
subplot(2,3,1)
plot((vert_means.SF_prod(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot((vert_means.MF_prod(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
title('Forage')
xlabel('Prod (g m^-^3 d^-^1)')

subplot(2,3,2)
plot((vert_means.SP_prod(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot((vert_means.MP_prod(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot((vert_means.LP_prod(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
title('Lg pel')
xlabel('Prod (g m^-^3 d^-^1)')

subplot(2,3,3)
plot((vert_means.SD_prod(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
title('Demersal')
xlabel('Prod (g m^-^3 d^-^1)')

subplot(2,3,4)
plot(log10(vert_means.SF_prod),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vert_means.MF_prod),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
legend('SF','MF')
%legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
xlabel('log_1_0 Prod (g m^-^3 d^-^1)')
ylabel('Depth (m)')

subplot(2,3,5)
plot(log10(vert_means.SP_prod),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vert_means.MP_prod),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vert_means.LP_prod),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
legend('SP','MP','LP')
legend('location','east')
xlabel('log_1_0 Prod (g m^-^3 d^-^1)')
ylabel('Depth (m)')

subplot(2,3,6)
plot(log10(vert_means.SD_prod),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
plot(log10(vert_means.MD_prod),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
plot(log10(vert_means.LD_prod),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
legend('SD','MD','LD')
legend('location','east')
xlabel('log_1_0 Prod (g m^-^3 d^-^1)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath mod '_prod_vert_mean_fntypes.png'])



%% Maps

%% All small
figure(9)
% SF
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.SF_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean prod SF (g m^-^2 d^-^1)')

% SD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.SD_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
set(gcf,'renderer','painters')
title('log10 mean prod SD (g m^-^2 d^-^1)')

% SP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.SP_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
set(gcf,'renderer','painters')
title('log10 mean prod SP (g m^-^2 d^-^1)')

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
print('-dpng',[ppath exper '_global_small_prod_subplot.png'])

%% All med
figure(10)
% MF
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MF_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean prod MF (g m^-^2 d^-^1)')

% MD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MD_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
set(gcf,'renderer','painters')
title('log10 mean prod MD (g m^-^2 d^-^1)')

% MP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.MP_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
set(gcf,'renderer','painters')
title('log10 mean prod MP (g m^-^2 d^-^1)')

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
print('-dpng',[ppath exper '_global_med_prod_subplot.png'])

%% All lg
figure(11)
% SF
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(LAT,LON,log10(spat_vert_ints.SF_prod))
% cmocean('matter')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-6 -3]);
% colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
% set(gcf,'renderer','painters')
% title('log10 mean All F (g m^-^2)')

% LD
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.LD_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
set(gcf,'renderer','painters')
title('log10 mean prod LD (g m^-^2 d^-^1)')

% LP
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(spat_vert_ints.LP_prod))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-6 -3]);
set(gcf,'renderer','painters')
title('log10 mean prod LP (g m^-^2 d^-^1)')

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
print('-dpng',[ppath exper '_global_lrg_prod_subplot.png'])

%% Flux times biomass
% %Rec