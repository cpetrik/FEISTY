% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files
% Map same quantities and scales as Fish-MIP protocol eg

clear all
close all
warning off

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

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

% plot info
[ni,nj]=size(lon);
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

cmS=cbrewer('div','Spectral',11);
cmS2=cmS(6:end,:);
cmS2=flipud(cmS2);

cmRB=cbrewer('div','RdYlBu',40);
cmRB=flipud(cmRB);
cmPO=cbrewer('div','PuOr',40);
cmPO=flipud(cmPO);

% save('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat',...
%     'mycmap','-append')
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmPR = mycmap;

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

%total
Zsf(ID)=sf_tot;
Zsp(ID)=sp_tot;
Zsd(ID)=sd_tot;
Zmf(ID)=mf_tot;
Zmp(ID)=mp_tot;
Zmd(ID)=md_tot;
Zlp(ID)=lp_tot;
Zld(ID)=ld_tot;
Zb(ID)=b_tot;

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

%sum over year to be per year instead of per day
Psf(ID)=sf_tprod;
Psp(ID)=sp_tprod;
Psd(ID)=sd_tprod;
Pmf(ID)=mf_tprod;
Pmp(ID)=mp_tprod;
Pmd(ID)=md_tprod;
Plp(ID)=lp_tprod;
Pld(ID)=ld_tprod;

ocean=NaN*ones(ni,nj);
ocean(ID)=ones(size(sf_tot));


%% Diff maps of all fish
%total consumer biomass kgC m-2

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
%from wet weight to g C
% All = All/9;
% PAll = PAll/9;

%% All Biom
figure(1)
subplot('Position',[0 0.6 0.5 0.32])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-180 180],'frame','on',...
    'Grid','off','FLineWidth',1)%,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(All))
colormap(cmPR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1e4]);
hcb = colorbar('v');
set(gcf,'renderer','painters')
title('Climatology Biomass All fishes (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_All_biomass_FishMIPscale.png'])

%% All production
figure(2)
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(PAll))
colormap(cmPR)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1e4]);
hcb = colorbar('v');
set(gcf,'renderer','painters')
title('Climatology Production All fishes (g m^-^2)')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_All_prod_FishMIPscale.png'])

%% All Biom
figure(5)
axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',[-180 180])
surfm(geolat_t,geolon_t,(All));
%xlim([-180 180])
%ylim([-90 90])
colormap(cmPR)
caxis([0 1e4]);
colorbar('v')
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%title('Climatology Biomass All fishes (g m^-^2)')
%stamp([harv '_' cfile])
print('-dpng',[ppath 'Climatol_' harv '_global_All_biomass_FishMIPscale_axes.png'])

