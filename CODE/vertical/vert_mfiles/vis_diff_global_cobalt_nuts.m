% Plot differences between COBALT only and COBALT-FEISTY 
% COBALT nuts, O2, chl vars

clear
close all

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



