% Save output of FEISTY Hindcast averaged over biomes
% for foodweb diagram
% 145 years, monthly 
% Use Zprod instead of Zloss
% Use prod instead of consumption (b/c don't have consuption)

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([bpath 'COBALT_biomes_50yr_means_5100.mat'],'biome_hist');

load([cpath 'hindcast_gridspec.mat'],'AREA_OCN','geolon_t','geolat_t');
AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);


%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
lpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

% Hindcast & Forecast together
load([fpath 'Means_hist_fore_',harv,'_cobalt_' cfile '.mat']);

zprod_hist = mzprod_hist + lzprod_hist;
zloss_hist = mzloss_hist + lzloss_hist;

%% Zoop, det, bent
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'mz_mean_hist','lz_mean_hist');
grid = csvread([cpath 'grid_csv.csv']);

% molN/m2 --> g/m2
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mz_mean_hist = mz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
lz_mean_hist = lz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
z_mean_hist = mz_mean_hist + lz_mean_hist;

%% Biomes
bios = nan(3,6);
flux = nan(3,7);
for B=1:3
    bid = find(biome_hist==B);
    bios(B,1) = nanmean(npp_hist(bid));
    bios(B,2) = nanmean(z_mean_hist(bid));
    bios(B,3) = nanmean(hF(bid));
    bios(B,4) = nanmean(hP(bid));
    bios(B,5) = nanmean(hB(bid));
    bios(B,6) = nanmean(hD(bid));
    
    flux(B,1) = nanmean(zprod_hist(bid));
    flux(B,2) = nanmean(zloss_hist(bid));
    flux(B,3) = nanmean(hFprod(bid));
    flux(B,4) = nanmean(det_hist(bid));
    flux(B,5) = nanmean(hDprod(bid));
    flux(B,6) = nanmean(hPprod(bid));
    flux(B,7) = nanmean(hFprod(bid));
end

%% save
bnames = {'LC','ECCS','ECSS'};
save([fpath 'Hist_' harv '_biomes_biom_flux_KKfoodweb_prods.mat'],'bnames',...
    'bios','flux');
save([lpath 'Hist_' harv '_biomes_biom_flux_KKfoodweb_prods.mat'],'bnames',...
    'bios','flux');




