% Visualize output of Spinup
% 50 years with movement, 10 years without dB/dt
% Compare 1st 10 years of movement & physiol and
% movement without dB/dt update

clear
close all

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90;
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5];...      %med grey


set(groot,'defaultAxesColorOrder',cm10);

%% with bio
exper = 'Spinup1988_move_ingest_v21_dt6h_All_fish03_';

load([fpath 'Means_' exper cfile '.mat']);

bF = sf_tmean+mf_tmean;
bP = sp_tmean+mp_tmean+lp_tmean;
bD = sd_tmean+md_tmean+ld_tmean;
bB = b_tmean;

bSFt = sf_tmean;
bMFt = mf_tmean;
bSPt = sp_tmean;
bMPt = mp_tmean;
bLPt = lp_tmean;
bSDt = sd_tmean;
bMDt = md_tmean;
bLDt = ld_tmean;
bBt  = b_tmean;

bSFm = mean(sf_smean(:,1:10),2);
bMFm = mean(mf_smean(:,1:10),2);
bSPm = mean(sp_smean(:,1:10),2);
bMPm = mean(mp_smean(:,1:10),2);
bLPm = mean(lp_smean(:,1:10),2);
bSDm = mean(sd_smean(:,1:10),2);
bMDm = mean(md_smean(:,1:10),2);
bLDm = mean(ld_smean(:,1:10),2);
bBm  = mean(b_smean(:,1:10),2);

clear sf_smean sp_smean sd_smean mf_smean mp_smean md_smean b_smean lp_smean ld_smean
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean b_tmean lp_tmean ld_tmean 
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean b_mean

%% no bio
exper2 = 'Spinup1988_move_prey_v23_nobio_dt6h_All_fish03_';

load([fpath 'Means_' exper2 cfile '.mat']);

nF = sf_tmean+mf_tmean;
nP = sp_tmean+mp_tmean+lp_tmean;
nD = sd_tmean+md_tmean+ld_tmean;
nB = b_tmean;

nSFt = sf_tmean;
nMFt = mf_tmean;
nSPt = sp_tmean;
nMPt = mp_tmean;
nLPt = lp_tmean;
nSDt = sd_tmean;
nMDt = md_tmean;
nLDt = ld_tmean;
nBt  = b_tmean;

nSFm = sf_mean;
nMFm = mf_mean;
nSPm = sp_mean;
nMPm = mp_mean;
nLPm = lp_mean;
nSDm = sd_mean;
nMDm = md_mean;
nLDm = ld_mean;
nBm  = b_mean;

clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean b_tmean lp_tmean ld_tmean 
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean b_mean

%% Does no dB/dt conserve mass?
nF(end)-nF(1)
nP(end)-nP(1)
nD(end)-nD(1)


%% Plots in time
y = time/12;

% All size classes of all together
figure(1)
plot(y,log10(bB(1:120)),'color',[0.75 0.75 0.75],'Linewidth',2); hold on;
plot(y,log10(bF(1:120)),'r','Linewidth',2); hold on;
plot(y,log10(bP(1:120)),'b','Linewidth',2); hold on;
plot(y,log10(bD(1:120)),'k','Linewidth',2); hold on;
plot(y,log10(nB),'color',[0.25 0.25 0.25],'Linewidth',2); hold on;
plot(y,log10(nF),'m','Linewidth',2); hold on;
plot(y,log10(nP),'c','Linewidth',2); hold on;
plot(y,log10(nD),'g','Linewidth',2); hold on;
legend('bB','bF','bP','bD','nB','nF','nP','nD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-2 2])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup ingest with & without dB/dt')
stamp(exper2)
print('-dpng',[ppath exper2 'all_types_comp.png'])

%% 
figure(2)
plot(y,log10(bBt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bSFt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bMFt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bSPt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bMPt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bLPt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bSDt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bMDt(1:120)),'Linewidth',1); hold on;
plot(y,log10(bLDt(1:120)),'Linewidth',1); hold on;
plot(y,log10(nBt),'--','Linewidth',1); hold on;
plot(y,log10(nSFt),'--','Linewidth',1); hold on;
plot(y,log10(nMFt),'--','Linewidth',1); hold on;
plot(y,log10(nSPt),'--','Linewidth',1); hold on;
plot(y,log10(nMPt),'--','Linewidth',1); hold on;
plot(y,log10(nLPt),'--','Linewidth',1); hold on;
plot(y,log10(nSDt),'--','Linewidth',1); hold on;
plot(y,log10(nMDt),'--','Linewidth',1); hold on;
plot(y,log10(nLDt),'--','Linewidth',1); hold on;
legend('bB','bSF','bMF','bSP','bMP','bLP','bSD','bMD','bLD','nB','nSF','nMF','nSP','nMP','nLP','nSD','nMD','nLD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-3 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup ingest with & without dB/dt')
stamp(exper2)
print('-dpng',[ppath exper2 'all_sizes.png'])

%% Diffs
figure(3)
plot(y,nB-bB(1:120),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,nF-bF(1:120),'r','Linewidth',2); hold on;
plot(y,nP-bF(1:120),'b','Linewidth',2); hold on;
plot(y,nD-bD(1:120),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-2 2])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2) difference without - with')
title('Spinup ingest with & without dB/dt')
stamp(exper2)
print('-dpng',[ppath exper2 'all_types_diff.png'])

%% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(GRD.ID)=bSFm;
Zsp(GRD.ID)=bSPm;
Zsd(GRD.ID)=bSDm;
Zmf(GRD.ID)=bMFm;
Zmp(GRD.ID)=bMPm;
Zmd(GRD.ID)=bMDm;
Zlp(GRD.ID)=bLPm;
Zld(GRD.ID)=bLDm;
Zb(GRD.ID)=bBm;

zAll = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
zAllF = Zsf+Zmf;
zAllP = Zsp+Zmp+Zlp;
zAllD = Zsd+Zmd+Zld;

Nsf=NaN*ones(ni,nj);
Nsp=NaN*ones(ni,nj);
Nsd=NaN*ones(ni,nj);
Nmf=NaN*ones(ni,nj);
Nmp=NaN*ones(ni,nj);
Nmd=NaN*ones(ni,nj);
Nlp=NaN*ones(ni,nj);
Nld=NaN*ones(ni,nj);
Nb=NaN*ones(ni,nj);

Nsf(GRD.ID)=nSFm;
Nsp(GRD.ID)=nSPm;
Nsd(GRD.ID)=nSDm;
Nmf(GRD.ID)=nMFm;
Nmp(GRD.ID)=nMPm;
Nmd(GRD.ID)=nMDm;
Nlp(GRD.ID)=nLPm;
Nld(GRD.ID)=nLDm;
Nb(GRD.ID)=nBm;

nAll = Nsp+Nsf+Nsd+Nmp+Nmf+Nmd+Nlp+Nld;
nAllF = Nsf+Nmf;
nAllP = Nsp+Nmp+Nlp;
nAllD = Nsd+Nmd+Nld;


%% Diff - bent
figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Nb - Zb))
cmocean('balance')
load coastlines;
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-25 25]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('benthic biomass (g m^-^2) diff without - with')
stamp(exper2)
print('-dpng',[ppath exper2 'global_BENT.png'])

%% Diffs
figure(5)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,nAllF-zAllF)
cmocean('balance')
load coastlines;
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('All F (g m^-^2) diff without - with')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,nAllD-zAllD)
cmocean('balance')
load coastlines;
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
set(gcf,'renderer','painters')
title('All D (g m^-^2) diff without - with')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,nAllP-zAllP)
cmocean('balance')
load coastlines;
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
set(gcf,'renderer','painters')
title('All P (g m^-^2) diff without - with')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,nAll-zAll)
cmocean('balance')
load coastlines;
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-10 10]);
set(gcf,'renderer','painters')
title('All fishes (g m^-^2) diff without - with')
stamp(exper2)
print('-dpng',[ppath exper2 'All_subplot.png'])




