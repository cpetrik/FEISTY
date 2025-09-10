% Visualize diff of movement with and without advection if nothing
% favorable outside current grid cell
% 1 year
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
exper = 'Spinup1988_move_enc_';
load([fpath exper 'v21_dt6h_All_fish03_Y1.mat']);

%%
%Time
sp_tmean=mean(S_Sml_p,1);
sf_tmean=mean(S_Sml_f,1);
sd_tmean=mean(S_Sml_d,1);
mp_tmean=mean(S_Med_p,1);
mf_tmean=mean(S_Med_f,1);
md_tmean=mean(S_Med_d,1);
lp_tmean=mean(S_Lrg_p,1);
ld_tmean=mean(S_Lrg_d,1);
b_tmean=mean(S_Bent_bio,1);

% Space
%whole yr
sp_mean=mean(S_Sml_p,2);
sf_mean=mean(S_Sml_f,2);
sd_mean=mean(S_Sml_d,2);
mp_mean=mean(S_Med_p,2);
mf_mean=mean(S_Med_f,2);
md_mean=mean(S_Med_d,2);
lp_mean=mean(S_Lrg_p,2);
ld_mean=mean(S_Lrg_d,2);
b_mean=mean(S_Bent_bio,2);

%%
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

Zsf(GRD.ID)=sf_mean;
Zsp(GRD.ID)=sp_mean;
Zsd(GRD.ID)=sd_mean;
Zmf(GRD.ID)=mf_mean;
Zmp(GRD.ID)=mp_mean;
Zmd(GRD.ID)=md_mean;
Zlp(GRD.ID)=lp_mean;
Zld(GRD.ID)=ld_mean;
Zb(GRD.ID)=b_mean;

sAll = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
sAllF = Zsf+Zmf;
sAllP = Zsp+Zmp+Zlp;
sAllD = Zsd+Zmd+Zld;

%% Not passive
exper = 'Spinup1988_move_enc_';
load([fpath exper 'v22_notpassive_dt6h_All_fish03_Y1.mat']);

%%
%Time
sp_tmean=mean(S_Sml_p,1);
sf_tmean=mean(S_Sml_f,1);
sd_tmean=mean(S_Sml_d,1);
mp_tmean=mean(S_Med_p,1);
mf_tmean=mean(S_Med_f,1);
md_tmean=mean(S_Med_d,1);
lp_tmean=mean(S_Lrg_p,1);
ld_tmean=mean(S_Lrg_d,1);
b_tmean=mean(S_Bent_bio,1);

% Space
%whole yr
sp_mean=mean(S_Sml_p,2);
sf_mean=mean(S_Sml_f,2);
sd_mean=mean(S_Sml_d,2);
mp_mean=mean(S_Med_p,2);
mf_mean=mean(S_Med_f,2);
md_mean=mean(S_Med_d,2);
lp_mean=mean(S_Lrg_p,2);
ld_mean=mean(S_Lrg_d,2);
b_mean=mean(S_Bent_bio,2);

%%
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

Msf(GRD.ID)=sf_mean;
Msp(GRD.ID)=sp_mean;
Msd(GRD.ID)=sd_mean;
Mmf(GRD.ID)=mf_mean;
Mmp(GRD.ID)=mp_mean;
Mmd(GRD.ID)=md_mean;
Mlp(GRD.ID)=lp_mean;
Mld(GRD.ID)=ld_mean;
Mb(GRD.ID)=b_mean;

mAll = Msp+Msf+Msd+Mmp+Mmf+Mmd+Mlp+Mld;
mAllF = Msf+Mmf;
mAllP = Msp+Mmp+Mlp;
mAllD = Msd+Mmd+Mld;

%% Diffs
dAll = mAll - sAll;
dAllF = mAllF - sAllF;
dAllP = mAllP - sAllP;
dAllD = mAllD - sAllD;
dAllB = Mb - Zb;

dts_B=(mB-sB);
dts_F=(mF-sF);
dts_P=(mP-sP);
dts_D=(mD-sD);

%% Plots in time
y = 1:365;

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
legend('pB','pF','pP','pD','npB','npF','npP','npD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.5 0.5])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['Y1 spin passive vs. not passive'])
print('-dpng',[ppath exper 'passive_notpassive_all_types_ts.png'])

figure(2)
plot(y,(dts_B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,(dts_F),'r','Linewidth',2); hold on;
plot(y,(dts_P),'b','Linewidth',2); hold on;
plot(y,(dts_D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-0.03 0.03])
xlabel('Time (y)')
ylabel('Biomass difference (g m^-^2)')
title(['Y1 spin Not passive - Passive'])
print('-dpng',[ppath exper 'passive_notpassive_diff_all_types_ts.png'])

%% bent
figure(3)
subplot('Position',[0 0.53 0.5 0.5])
%No move
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log_1_0 Passive Benthos')

%Not passive
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Mb))
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 2]);
set(gcf,'renderer','painters')
title('log_1_0 Not passive Benthos')

%Diff
subplot('Position',[0.25 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,dAllB)
cmocean('balance')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.1 0.1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Not passive - Passive')
print('-dpng',[ppath exper 'passive_notpassive_diff_BENT.png'])

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
title('All F (g m^-^2) Not passive - Passive')

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
title('All D (g m^-^2)')

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
title('All P (g m^-^2)')

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
title('All fishes (g m^-^2)')
stamp('')
print('-dpng',[ppath exper 'passive_notpassive_diff_All_subplot.png'])

