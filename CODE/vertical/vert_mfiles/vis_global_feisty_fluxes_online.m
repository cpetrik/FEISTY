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
load([fpath '19900101.ocean_feisty_fluxes_z_means.mat'])

%%
tF = tSF+tMF;
tP = tSP+tMP+tLP;
tD = tSD+tMD+tLD;

sF = sSF+sMF;
sP = sSP+sMP+sLP;
sD = sSD+sMD+sLD;
sAll = sF+sP+sD;

FracPD = sP ./ (sP + sD);
FracPF = sP ./ (sP + sF);
FracLM = (sLP+sLD) ./ (sLP+sLD+sMF+sMP+sMD);

vF = vSF+vMF;
vP = vSP+vMP+vLP;
vD = vSD+vMD+vLD;
vAll = vF+vP+vD;

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
y = 1:12; %120;

% All size classes of all
figure(1)
plot(y,log10(tSF(y)),'Linewidth',1); hold on;
plot(y,log10(tMF(y)),'Linewidth',1); hold on;
plot(y,log10(tSP(y)),'Linewidth',1); hold on;
plot(y,log10(tMP(y)),'Linewidth',1); hold on;
plot(y,log10(tLP(y)),'Linewidth',1); hold on;
plot(y,log10(tSD(y)),'Linewidth',1); hold on;
plot(y,log10(tMD(y)),'Linewidth',1); hold on;
plot(y,log10(tLD(y)),'Linewidth',1); hold on;
plot(y,log10(tBE(y)),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD','B')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Integrated Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath mod '_ts_logmean_feisty_all_sizes.png'])

% Fn Types
figure(2)
plot(y,log10(tBE(y)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(tF(y)),'r','Linewidth',2); hold on;
plot(y,log10(tP(y)),'b','Linewidth',2); hold on;
plot(y,log10(tD(y)),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log10 Integrated Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath mod '_ts_logmean_feisty_all_types.png'])


figure(3)
subplot(3,3,1)
plot(y,(tSF(y)),'color',cm10(1,:),'Linewidth',1); 
subplot(3,3,4)
plot(y,(tMF(y)),'color',cm10(2,:),'Linewidth',1); 
subplot(3,3,2)
plot(y,(tSP(y)),'color',cm10(3,:),'Linewidth',1); 
subplot(3,3,5)
plot(y,(tMP(y)),'color',cm10(4,:),'Linewidth',1); 
subplot(3,3,8)
plot(y,(tLP(y)),'color',cm10(5,:),'Linewidth',1); 
subplot(3,3,3)
plot(y,(tSD(y)),'color',cm10(6,:),'Linewidth',1); 
subplot(3,3,6)
plot(y,(tMD(y)),'color',cm10(7,:),'Linewidth',1); 
subplot(3,3,9)
plot(y,(tLD(y)),'color',cm10(8,:),'Linewidth',1);
subplot(3,3,7)
plot(y,(tBE(y)),'color',cm10(9,:),'Linewidth',1);
%xlim([y(1) y(end)])
xlabel('Time (mo)')
ylabel('Integrated Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath mod '_ts_mean_feisty_all_sizes.png'])

% Fn Types
figure(4)
subplot(2,2,1)
plot(y,(tF(y)),'r','Linewidth',2);
subplot(2,2,2)
plot(y,(tP(y)),'b','Linewidth',2); 
subplot(2,2,3)
plot(y,(tD(y)),'k','Linewidth',2); hold on;
subplot(2,2,4)
plot(y,(tBE(y)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
%xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Integrated Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath mod '_ts_mean_feisty_all_types.png'])

%% Vert distrib
figure(5)
subplot(1,2,1)
plot(log10(vSF(:,1)),-1*z_l,'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vMF(:,1)),-1*z_l,'color',cm10(2,:),'Linewidth',1); hold on;
plot(log10(vSP(:,1)),-1*z_l,'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vMP(:,1)),-1*z_l,'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vLP(:,1)),-1*z_l,'color',cm10(5,:),'Linewidth',1); hold on;
plot(log10(vSD(:,1)),-1*z_l,'color',cm10(6,:),'Linewidth',1); hold on;
% plot(log10(vMD(:,1)),-1*z_l,'color',cm10(7,:),'Linewidth',1); hold on;
% plot(log10(vLD(:,1)),-1*z_l,'color',cm10(8,:),'Linewidth',1); hold on;
% plot(log10(vBE(:,1)),-1*z_l,'color',cm10(9,:),'Linewidth',1); hold on;
% legend('SF','MF','SP','MP','LP','SD','MD','LD','B')
legend('SF','MF','SP','MP','LP','SD')
legend('location','east')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')

subplot(1,2,2)
plot(log10(vF(:,1)),-1*z_l,'r','Linewidth',1); hold on;
plot(log10(vP(:,1)),-1*z_l,'b','Linewidth',1); hold on;
plot(log10(vD(:,1)),-1*z_l,'k','Linewidth',1); hold on;
legend('F','P','D')
legend('location','east')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath exper '_vert_mean_feisty_subplot.png'])

%% Vert distrib - upper ocean
figure(6)
subplot(1,2,1)
plot(log10(vSF(1:10,1)),-1*z_l(1:10),'color',cm10(1,:),'Linewidth',1); hold on;
plot(log10(vMF(1:10,1)),-1*z_l(1:10),'color',cm10(2,:),'Linewidth',1); hold on;
plot(log10(vSP(1:10,1)),-1*z_l(1:10),'color',cm10(3,:),'Linewidth',1); hold on;
plot(log10(vMP(1:10,1)),-1*z_l(1:10),'color',cm10(4,:),'Linewidth',1); hold on;
plot(log10(vLP(1:10,1)),-1*z_l(1:10),'color',cm10(5,:),'Linewidth',1); hold on;
plot(log10(vSD(1:10,1)),-1*z_l(1:10),'color',cm10(6,:),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD')
legend('location','southeast')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')

subplot(1,2,2)
plot(log10(vF(1:10,1)),-1*z_l(1:10),'r','Linewidth',1); hold on;
plot(log10(vP(1:10,1)),-1*z_l(1:10),'b','Linewidth',1); hold on;
plot(log10(vD(1:10,1)),-1*z_l(1:10),'k','Linewidth',1); hold on;
legend('F','P','D')
legend('location','southeast')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath exper '_vert_upper_mean_feisty_subplot.png'])

%% Maps
% bent
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sBE))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean benthic biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_Bent.png'])

%% All 4 on subplots
figure(8)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sF))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sD))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sP))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sAll))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
%stamp(cfile)
print('-dpng',[ppath exper '_global_All_subplot.png'])

%% All 4 on subplots
figure(18)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sF))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sD))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sP))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,log10(sAll))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
%stamp(cfile)
print('-dpng',[ppath exper '_global_All_subplot_newcb.png'])

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(9)
subplot('Position',[0 0.53 0.5 0.5])
%P:D
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,FracPD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%P:F
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,FracPF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%L:M
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(LAT,LON,FracLM)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Fraction Large vs. Medium')
stamp(cfile)
print('-dpng',[ppath exper '_global_ratios_subplot.png'])
