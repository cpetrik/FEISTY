% Visualize output of FEISTY with Arctic projection
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);


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

Zsf(ID)=sf_mean;
Zsp(ID)=sp_mean;
Zsd(ID)=sd_mean;
Zmf(ID)=mf_mean;
Zmp(ID)=mp_mean;
Zmd(ID)=md_mean;
Zlp(ID)=lp_mean;
Zld(ID)=ld_mean;
Zb(ID)=b_mean;

% Diff groups
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;

MFP = Zmf+Zmp;
LLP = Zlp;

%% Save for Jerome
fish_biomass_units = 'g wet weight m-2';
length_units = 'mm';
weight_units = 'g wet weight';
size_bins = {'min','geometric mean','max'};
S_lengths = [0.001,0.02,0.5];
M_lengths = [0.5,11.2,250];
L_lengths = [250,5600,125000];
S_weights = [4.6,13,36.8];
M_weights = [36.8,104,292.4];
L_weights = [292.4,824,2320.8];

save([fpath 'FEISTY_biomass_means_1deg_ESM26_5yr_clim_191_195.mat'],...
    'All','AllF','AllP','AllD','AllS','AllM','AllL','fish_biomass_units',...
    'length_units','weight_units','size_bins','S_lengths','M_lengths',...
    'L_lengths','S_weights','M_weights','L_weights','geolat_t','geolon_t',...
    'MFP','LLP')

