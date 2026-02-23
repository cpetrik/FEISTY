% Plot percent differences between COBALT only and COBALT-FEISTY 
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
NtNH4 = tNH4(1:60);
NtNO3 = tNO3(1:60);
NtO2  = tO2(1:60);
NtCHL = tCHL(1:60);

NvNH4 = vNH4(:,1);
NvNO3 = vNO3(:,1);
NvO2  = vO2(:,1);
NvCHL = vCHL(:,1);

NsNH4 = sNH4(:,:,1);
NsNO3 = sNO3(:,:,1);
NsO2 = sO2(:,:,1);
NsCHL = sCHL(:,:,1);
NsCHLs = sCHLs(:,:,1);

clear tNH4 tNO3 vNH4 vNO3 sNH4 sNO3
clear tO2 tCHL vO2 vCHL sO2 sCHL sCHLs

%% OFFLINE -----------------------------------------------------------
%fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
fpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/COBALTonly/';

mod = 'OM4_05_COBALTv3_FEISTYoff';

load([fpath 'ocean_cobalt_nuts_month_z.199001-199412_means.mat'])

%% 
FtNH4 = tNH4(1:60);
FtNO3 = tNO3(1:60);
FtO2  = tO2(1:60);
FtCHL = tCHL(1:60);

FvNH4 = vNH4(:,1);
FvNO3 = vNO3(:,1);
FvO2  = vO2(:,1);
FvCHL = vCHL(:,1);

FsNH4 = sNH4(:,:,1);
FsNO3 = sNO3(:,:,1);
FsO2  = sO2(:,:,1);
FsCHL = sCHL(:,:,1);
FsCHLs = sCHLs(:,:,1);

clear tNH4 tNO3 vNH4 vNO3 sNH4 sNO3
clear tO2 tCHL vO2 vCHL sO2 sCHL sCHLs

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

dtNH4 = 100*(NtNH4 - FtNH4) ./ FtNH4;
dtNO3 = 100*(NtNO3 - FtNO3) ./ FtNO3;
dtO2 = 100*(NtO2 - FtO2) ./ FtO2;
dtCHL = 100*(NtCHL - FtCHL) ./ FtCHL;

dsNH4 = 100*(NsNH4 - FsNH4) ./ FsNH4;
dsNO3 = 100*(NsNO3 - FsNO3) ./ FsNO3;
dsO2 = 100*(NsO2 - FsO2) ./ FsO2;
dsCHL = 100*(NsCHL - FsCHL) ./ FsCHL;
dsCHLs = 100*(NsCHLs - FsCHLs) ./ FsCHLs;

dvNH4 = 100*(NvNH4 - FvNH4) ./ FvNH4;
dvNO3 = 100*(NvNO3 - FvNO3) ./ FvNO3;
dvO2 = 100*(NvO2 - FvO2) ./ FvO2;
dvCHL = 100*(NvCHL - FvCHL) ./ FvCHL;

tmos = 1:60;

%% PLOTS OF PercentDIFFS ---------------------------------------
% Time series
% NO Log10
figure(7)
subplot(3,1,1)
plot(tmos,(dtNH4(tmos)+eps),'color',cm10(9,:),'LineWidth',2); hold on;
plot(tmos,(dtNO3(tmos)+eps),'color',cm10(10,:),'LineWidth',2); hold on;
title('Integrated Concentration (mol m^-^2) %Difference Online - Offline')
legend({'NH4','NO3'})
legend('location','northwest')
subplot(3,1,2)
plot(tmos,(dtCHL(tmos)+eps),'color',cm10(2,:),'LineWidth',2); hold on;
ylabel('Online - Offline')
title('Integrated Chl Concentration (ug m^-^2) %Difference Online - Offline')
subplot(3,1,3)
plot(tmos,(dtO2(tmos)+eps),'color',cm10(1,:),'LineWidth',2); 
title('O2 Concentration (mol m^-^2) %Difference Online - Offline')
xlabel('Months')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_ts_Pdiff_cobalt_nuts.png'])

%% Vert
figure(8)
subplot(1,2,1)
plot((dvNH4(1:9,1)),-1*z_l(1:9),'color',cm10(9,:),'LineWidth',2); hold on; 
plot((dvNO3(1:9,1)),-1*z_l(1:9),'color',cm10(10,:),'LineWidth',2); hold on;
plot((dvO2(1:9,1)),-1*z_l(1:9),'color',cm10(1,:),'LineWidth',2); hold on;
legend({'NH4','NO3','O2'})
legend('location','southeast')
title('Mean Concen (mol m^-^3) %Difference Online - Offline')
subplot(1,2,2)
plot((dvCHL(1:9,1)),-1*z_l(1:9),'color',cm10(2,:),'LineWidth',2); hold on;
title('Mean Chl Concen (ug m^-^3) %Difference Online - Offline')
ylabel('Depth (m)')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_vert_Pdiff_cobalt_nuts.png'])


%% Maps
% NO3, NH4
figure(9)
subplot(1,2,1) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNO3(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NO3 %difference (mol m^2)')
text(4,2.5,'Online - Offline','HorizontalAlignment','center')

subplot(1,2,2) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNH4(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NH4 %difference (mol m^2)')

print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_no3_nh4.png'])

%% CHLs
figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsCHLs(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Surf Chl %difference (m^3) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_Pdiff_surfChl.png'])

%% Chl
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsCHL(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('Int Chl %difference (m^2) Online - Offline')
stamp('')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_mean_Pdiff_intChl.png'])

%% O2
figure(12) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsO2(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('O2 %difference (mol m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_O2.png'])

%% NH4
figure(13)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNH4(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NH4 %difference (mol m^2) Online - Offline')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_NH4.png'])

%% NO3
figure(14)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,(squeeze(dsNO3(:,:,1))))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-100 100]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('NO3 %difference (mol m^2)')
print('-dpng',[ppath 'OM4_05_COBALTv3_FEISTY_on-off_map_Pdiff_NO3.png'])



