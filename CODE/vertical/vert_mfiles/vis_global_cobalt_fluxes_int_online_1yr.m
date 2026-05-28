% Look at COBALT vint fluxes from online sim
% 12 mo of output

clear
close all

%%
fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

cfile ='NoDc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODe/Figs/Matlab_New_sizes/',...
    cfile,'/Cobalt_Feisty/'];

exper = 'OM4_05_COBALTv3_FEISTYon_021326';
Yr = '1990';
mod = [exper '_' Yr];

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
load([fpath Yr '0101.ocean_cobalt_fluxes_int_means.mat'])

%% molN/m2/s to gC/m2/d
S2D= 60*60*24;
NtoC= (106/16) * 12.01;

N2CD = NtoC * S2D;

% jingest_n_hp_100
% Size:       720x576x12
% missing_value = 1.000000020040877e+20
% units         = 'mol m-2 s-1'
% long_name     = 'Higher predator ingestion of nitrogen integral in upper 100m'


% jhploss_nlgz_100
% Size:       720x576x12
% missing_value = 1.000000020040877e+20
% units         = 'mol m-2 s-1'
% long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m'

%%
% z_l_ts = repmat(z_l,1,length(tmos));
% [z_l2,tts] = meshgrid(tmos,-1*z_l);

%% Time series
tmos = 1:12;

figure(1)
plot(tmos,log10(N2CD*tMZzl(tmos)+eps),'color',cm10(4,:)); hold on; 
plot(tmos,log10(N2CD*tHPinj(tmos)+eps),'color',cm10(5,:)); hold on; 
plot(tmos,log10(N2CD*tMZloss(tmos)+eps),'color',cm10(6,:)); hold on; 
plot(tmos,log10(N2CD*tLZloss(tmos)+eps),'color',cm10(7,:)); hold on;
legend({'MZzloss','HPingest','MZhploss','LZhploss'})
legend('location','eastoutside')
title('log_1_0 Integrated Flux (gC m^-^2 d^-^1)')
xlabel('Months')
stamp(Yr)
print('-dpng',[ppath mod '_ts_mean_cobalt_fluxes_int.png'])

figure(2)
plot(tmos,log10(N2CD*tHPinj+eps),'color',cm10(6,:)); hold on; 
plot(tmos,log10(N2CD*(tMZloss+tLZloss)+eps),'--k'); hold on; 
legend({'HPingest','MZ+LZhploss'})
legend('location','southeast')
title('log_1_0 Integrated Flux (gC m^-^2 d^-^1)')
xlabel('Months')
stamp(Yr)
print('-dpng',[ppath mod '_ts_mean_cobalt_hploss_ingest_int.png'])

%% Maps
% All 4 on subplots
figure(8)
% all F

% MZ zloss
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CD*sMZzl)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 MZ zloss vint (gC m^-^2 d^-^1)')

% MZ HPloss
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CD*sMZloss)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 MZ HPloss vint (gC m^-^2 d^-^1)')

% LZ HPloss
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CD*sLZloss)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 LZ HPloss vint (gC m^-^2 d^-^1)')

% HP ingest 
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,log10(squeeze(N2CD*sHPinj)))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 -1]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 HPingest vint (gC m^-^2 d^-^1)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Zloss_int_subplot.png'])

%% Adds up?
sdiff = sHPinj-sLZloss-sMZloss;
% HP ingest 
figure(7)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,sdiff)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('HPingest vint - MZloss - LZloss (gC m^-^2 d^-^1)')
stamp(Yr)
print('-dpng',[ppath mod '_global_Zloss_int_diff.png'])

%% HPloss vs Z biomass -> functional response
load([fpath 'ocean_cobalt_feisty_forcing_',Yr,'_means.mat'],'sTp','tTp',...
    'tMz','tLz','sMz','sLz')

sTdep = exp(0.063 * (sTp-10.0));
tTdep = exp(0.063 * (tTp-10.0));

sThpM = sMZloss ./ sTdep;
sThpL = sLZloss ./ sTdep;

tThpM = tMZloss ./ tTdep;
tThpL = tLZloss ./ tTdep;

%%
figure(9)
plot(tMz,tThpM,'.k','MarkerSize',15)
xlabel('MZ biomass')
ylabel('MZ HPloss normalized by temp-dep')
title('Global monthly means')
print('-dpng',[ppath exper Yr '_global_vint_MZvsHP_tmean'])

%%
figure(10)
plot(tLz,tThpL,'.k','MarkerSize',15)
xlabel('LZ biomass')
ylabel('LZ HPloss normalized by temp-dep')
title('Global monthly means')
print('-dpng',[ppath exper Yr '_global_vint_LZvsHP_tmean'])

%%
figure(11)
plot(sMz(:),sThpM(:),'.k','MarkerSize',10)
xlabel('MZ biomass')
ylabel('MZ HPloss normalized by temp-dep')
title('Vertically integrated, annual means')
print('-dpng',[ppath exper Yr '_global_vint_MZvsHP_smean.png'])

%%
figure(12)
plot(sLz(:),sThpL(:),'.k','MarkerSize',15)
xlabel('LZ biomass')
ylabel('LZ HPloss normalized by temp-dep')
title('Vertically integrated, annual means')
print('-dpng',[ppath exper Yr '_global_vint_LZvsHP_smean.png'])



