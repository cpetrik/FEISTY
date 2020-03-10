% Historic time period (1861-2005) at all locations
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_Historic_',harv,'_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% plot info
[ni,nj]=size(geolon_t);

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

Zsf(ID)=sf_mean5;
Zsp(ID)=sp_mean5;
Zsd(ID)=sd_mean5;
Zmf(ID)=mf_mean5;
Zmp(ID)=mp_mean5;
Zmd(ID)=md_mean5;
Zlp(ID)=lp_mean5;
Zld(ID)=ld_mean5;
Zb(ID)=b_mean5;

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

save([fpath 'FEISTY_biomass_means_ESM2M_hist9095.mat'],...
    'All','AllF','AllP','AllD','AllS','AllM','AllL','fish_biomass_units',...
    'length_units','weight_units','size_bins','S_lengths','M_lengths',...
    'L_lengths','S_weights','M_weights','L_weights','geolat_t','geolon_t',...
    'MFP','LLP')

