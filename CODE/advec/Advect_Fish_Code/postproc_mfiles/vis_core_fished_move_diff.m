% Visualize output of Spinup 
% 150 years
% Saved as mat files

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

%%
exper = 'CORE_Hindcast1988_no_move_';
load([fpath 'Means_' exper cfile '.mat']);

sF = sf_tmean+mf_tmean;
sP = sp_tmean+mp_tmean+lp_tmean;
sD = sd_tmean+md_tmean+ld_tmean;
sB = b_tmean;

% Plots in space
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(GRD.ID)=sf_mean20;
Zsp(GRD.ID)=sp_mean20;
Zsd(GRD.ID)=sd_mean20;
Zmf(GRD.ID)=mf_mean20;
Zmp(GRD.ID)=mp_mean20;
Zmd(GRD.ID)=md_mean20;
Zlp(GRD.ID)=lp_mean20;
Zld(GRD.ID)=ld_mean20;
Zb(GRD.ID)=b_mean20;

sAll = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
sAllF = Zsf+Zmf;
sAllP = Zsp+Zmp+Zlp;
sAllD = Zsd+Zmd+Zld;

%%
exper = 'CORE_Hindcast_move_prey_v5_';
load([fpath 'Means_' exper cfile '.mat']);

mF = sf_tmean+mf_tmean;
mP = sp_tmean+mp_tmean+lp_tmean;
mD = sd_tmean+md_tmean+ld_tmean;
mB = b_tmean;

Msf=NaN*ones(ni,nj);
Msp=NaN*ones(ni,nj);
Msd=NaN*ones(ni,nj);
Mmf=NaN*ones(ni,nj);
Mmp=NaN*ones(ni,nj);
Mmd=NaN*ones(ni,nj);
Mlp=NaN*ones(ni,nj);
Mld=NaN*ones(ni,nj);
Mb=NaN*ones(ni,nj);

Msf(GRD.ID)=sf_mean20;
Msp(GRD.ID)=sp_mean20;
Msd(GRD.ID)=sd_mean20;
Mmf(GRD.ID)=mf_mean20;
Mmp(GRD.ID)=mp_mean20;
Mmd(GRD.ID)=md_mean20;
Mlp(GRD.ID)=lp_mean20;
Mld(GRD.ID)=ld_mean20;
Mb(GRD.ID)=b_mean20;

mAll = Msp+Msf+Msd+Mmp+Mmf+Mmd+Mlp+Mld;
mAllF = Msf+Mmf;
mAllP = Msp+Mmp+Mlp;
mAllD = Msd+Mmd+Mld;

%% Diffs
dAll = sAll - mAll;
dAllF = sAllF - mAllF;
dAllP = sAllP - mAllP;
dAllD = sAllD - mAllD;
dAllB = Zb - Mb;

dts_B=(sB-mB);
dts_F=(sF-mF);
dts_P=(sP-mP);
dts_D=(sD-mD);

%% Plots in time
y = (time/12)+1988;

% All size classes of all
figure(1)
plot(y,log10(sB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(sF),'r','Linewidth',2); hold on;
plot(y,log10(sP),'b','Linewidth',2); hold on;
plot(y,log10(sD),'k','Linewidth',2); hold on;
plot(y,log10(mB),'--','color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(mF),'--r','Linewidth',2); hold on;
plot(y,log10(mP),'--b','Linewidth',2); hold on;
plot(y,log10(mD),'--k','Linewidth',2); hold on;
legend('sB','sF','sP','sD','mB','mF','mP','mD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.5 0.5])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['CORE stationary vs. move prey'])
%print('-dpng',[ppath 'all_types.png'])

figure(2)
plot(y,(sB-mB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(sF-mF),'r','Linewidth',2); hold on;
plot(y,(sP-mP),'b','Linewidth',2); hold on;
plot(y,(sD-mD),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.01 0.01])
xlabel('Time (y)')
ylabel('Biomass difference (g m^-^2)')
title(['CORE stationary vs. move prey'])

%% bent
figure(3)
subplot('Position',[0 0.53 0.5 0.5])
%No move
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Zb)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Demersals')

%Move
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,Mb)
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('Fraction Large Pelagics vs. Forage Fishes')

%Diff
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dAllB)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.5 0.5]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Stationary - Move prey')
%print('-dpng',[ppath exper 'global_BENT.png'])

%% ALL
figure(4)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dAllF)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('log10 mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dAllD)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
set(gcf,'renderer','painters')
title('log10 mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dAllP)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
set(gcf,'renderer','painters')
title('log10 mean All P (g m^-^2)')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dAll)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
set(gcf,'renderer','painters')
title('log10 mean All fishes (g m^-^2)')
stamp('')
%print('-dpng',[ppath exper 'global_All_subplot.png'])

