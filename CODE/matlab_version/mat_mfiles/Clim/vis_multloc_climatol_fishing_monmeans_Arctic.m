% Visualize output of FEISTY with Arctic projection
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

close all

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

geolat_t=lat;
geolon_t=lon;

% colors
cmBP=cbrewer('seq','BuPu',50);
cmGy=cbrewer('seq','Greys',50);
cmB=cbrewer('seq','Blues',50);

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

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;

mf_my(mf_my<0)=0;
Cmf(ID)=mf_my;
Cmp(ID)=mp_my;
Cmd(ID)=md_my;
Clp(ID)=lp_my;
Cld(ID)=ld_my;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_mean));

% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllM+AllL);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
FracPDs = Zsp ./ (Zsp+Zsd);
FracPDm = Zmp ./ (Zmp+Zmd);
FracPDl = Zlp ./ (Zlp+Zld);
FracPFs = Zsp ./ (Zsp+Zsf);
FracPFm = Zmp ./ (Zmp+Zmf);
FracPFvDs = (Zsp+Zsf) ./ (Zsp+Zsf+Zsd);
FracPFvDm = (Zmp+Zmf) ./ (Zmp+Zmf+Zmd);

%%
figure(1)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,log10(AllF));
shading interp
colorbar
colormap(cmB)
caxis([-2 2]);
m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
title('log_1_0 mean All F (g m^-^2)')
print('-dpng',[ppath 'Climatol_' harv '_arctic_F.png'])


%%
figure(2)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,log10(AllP));
shading interp
colorbar
colormap(cmB)
caxis([-2 2]);
m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
title('log_1_0 mean All P (g m^-^2)')
print('-dpng',[ppath 'Climatol_' harv '_arctic_P.png'])

%%
figure(3)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,log10(AllD));
shading interp
colorbar
colormap(cmB)
caxis([-2 2]);
m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
title('log_1_0 mean All D (g m^-^2)')
print('-dpng',[ppath 'Climatol_' harv '_arctic_D.png'])

%%
figure(4)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,log10(All));
shading interp
colorbar
colormap(cmB)
caxis([-2 2]);
m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
title('log_1_0 mean All Fishes (g m^-^2)')
print('-dpng',[ppath 'Climatol_' harv '_arctic_All.png'])
