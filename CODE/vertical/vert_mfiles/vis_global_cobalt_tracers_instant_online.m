% Look at COBALT-FEISTY tracers instant (nuts) from online sim

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';

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
load([fpath 'ocean_cobalt_tracers_instant.199001-199912_means.mat'])

%% flux in sec
S2D = 60*60*24;
D2Y = S2D*365;
N2C= (106/16) * 12.01;

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
tmos = 1:60;

figure(1)
subplot(3,1,1)
plot(tmos,log10(tDIC(tmos)+eps),'color',cm10(8,:)); hold on;
plot(tmos,log10(tDOC(tmos)+eps),'color',cm10(1,:)); hold on;
plot(tmos,log10(tPOC(tmos)+eps),'color',cm10(3,:)); hold on;
legend({'DIC','DOC','POC'})
legend('location','eastoutside')
ylabel('log_1_0 molC m^-^2')
title('Carbon')

subplot(3,1,2)
plot(tmos,log10(tO2(tmos)+eps),'color',cm10(5,:)); hold on;
ylabel('log_1_0 molO_2 m^-^2')
title('Oxygen')

subplot(3,1,2)
plot(tmos,log10(D2Y*N2C*tFN(tmos)+eps),'color',cm10(2,:)); hold on;
title('Total Sinking Flux')
title('log_1_0 gC m^-^2 y^-^1')
xlabel('Months')
stamp('')
print('-dpng',[ppath exper '_ts_mean_ocean_nuts_inst.png'])


%% Maps
% POC
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(12.01*sPOC(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-3 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 integrated POC (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_POCint.png'])

%% DOC
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(12.01*sDOC(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-3 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 integrated DOC (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_DOCint.png'])

%% DIC
figure(5)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(12.01*sDIC(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-1 0]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 integrated DIC  (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_DICint.png'])

%% O2
figure(6)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(sO2(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-1 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 integrated O_2 (mol m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_O2int.png'])

%% FnTot_100
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(D2Y*N2C*sFN(:,:,1))))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%clim([-0.5 0.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 total carbon flux (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath exper '_global_FnTot100.png'])

