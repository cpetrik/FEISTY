% Look at COBALT-FEISTY tracers (fish) from online sim

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

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
geolon = double(geolon);
geolat = double(geolat);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;  

%% 
load([fpath 'ocean_feisty_tracers_z.199001-199912_means.mat'])

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

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
y = 1:120;

figure(1)
% All size classes of all
figure(1)
plot(y,log10(tBE(y)),'Linewidth',1); hold on;
plot(y,log10(tSF(y)),'Linewidth',1); hold on;
plot(y,log10(tMF(y)),'Linewidth',1); hold on;
plot(y,log10(tSP(y)),'Linewidth',1); hold on;
plot(y,log10(tMP(y)),'Linewidth',1); hold on;
plot(y,log10(tLP(y)),'Linewidth',1); hold on;
plot(y,log10(tSD(y)),'Linewidth',1); hold on;
plot(y,log10(tMD(y)),'Linewidth',1); hold on;
plot(y,log10(tLD(y)),'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Integrated Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath mod '_ts_mean_feisty_all_sizes.png'])

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
print('-dpng',[ppath mod '_ts_mean_feisty_all_types.png'])


%% Vert distrib
figure(2)
subplot(1,2,1)
plot(log10(tBE(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tSF(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tMF(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tSP(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tMP(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tLP(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tSD(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tMD(:,1)),-1*z_l,'Linewidth',1); hold on;
plot(log10(tLD(:,1)),-1*z_l,'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','east')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')

subplot(1,2,2)
plot(log10(tBE(:,1)),-1*z_l,'Linewidth',1,'color',[0.5 0.5 0.5]); hold on;
plot(log10(tF(:,1)),-1*z_l,'Linewidth',1,'r'); hold on;
plot(log10(tP(:,1)),-1*z_l,'Linewidth',1,'b'); hold on;
plot(log10(tD(:,1)),-1*z_l,'Linewidth',1,'k'); hold on;
legend('B','F','P','D')
legend('location','east')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath exper '_vert_mean_feisty_subplot.png'])

%% Vert distrib - upper ocean
figure(10)
subplot(1,2,1)
plot(log10(tBE(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tSF(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tMF(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tSP(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tMP(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tLP(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tSD(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tMD(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
plot(log10(tLD(1:10,1)),-1*z_l(1:10),'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','east')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')

subplot(1,2,2)
plot(log10(tBE(1:10,1)),-1*z_l(1:10),'Linewidth',1,'color',[0.5 0.5 0.5]); hold on;
plot(log10(tF(1:10,1)),-1*z_l(1:10),'Linewidth',1,'r'); hold on;
plot(log10(tP(1:10,1)),-1*z_l(1:10),'Linewidth',1,'b'); hold on;
plot(log10(tD(1:10,1)),-1*z_l(1:10),'Linewidth',1,'k'); hold on;
legend('B','F','P','D')
legend('location','east')
title('log_1_0 Mean Biomass (g m^-^3)')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath exper '_vert_upper_mean_feisty_subplot.png'])

%% Maps
% bent
figure(3)
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
figure(4)
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

%% Ratios on subplots red-white-blue
% 3 figure subplot P:D, P:F, M:L
figure(5)
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
