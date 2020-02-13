% Save output of FEISTY Forecast averaged over biomes
% for foodweb diagram
% 95 years, monthly 
% Use Zprod instead of Zloss
% Use prod instead of consumption (b/c don't have consuption)

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([bpath 'COBALT_biomes_50yr_means_5100.mat'],'biome_fore');

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
zprod_fore = mzprod_fore + lzprod_fore;
zloss_fore = mzloss_fore + lzloss_fore;

%% Zoop, det, bent
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'mz_mean_fore','lz_mean_fore');
grid = csvread([cpath 'grid_csv.csv']);


% molN/m2 --> g/m2
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mz_mean_fore = mz_mean_fore * (106.0/16.0) * 12.01 * 9.0;
lz_mean_fore = lz_mean_fore * (106.0/16.0) * 12.01 * 9.0;
z_mean_fore = mz_mean_fore + lz_mean_fore;

%% Biomes
bios = nan(3,6);
flux = nan(3,7);
for B=1:3
    bid = find(biome_fore==B);
    bios(B,1) = nanmean(npp_fore(bid));
    bios(B,2) = nanmean(z_mean_fore(bid));
    bios(B,3) = nanmean(cF(bid));
    bios(B,4) = nanmean(cP(bid));
    bios(B,5) = nanmean(cB(bid));
    bios(B,6) = nanmean(cD(bid));
    
    flux(B,1) = nanmean(zprod_fore(bid));
    flux(B,2) = nanmean(zloss_fore(bid));
    flux(B,3) = nanmean(cFprod(bid));
    flux(B,4) = nanmean(det_fore(bid));
    flux(B,5) = nanmean(cDprod(bid));
    flux(B,6) = nanmean(cPprod(bid));
    flux(B,7) = nanmean(cFprod(bid));
end

%% save
bnames = {'LC','ECCS','ECSS'};
save([fpath 'Fore_' harv '_biomes_biom_flux_KKfoodweb_prods.mat'],'bnames',...
    'bios','flux');
save([lpath 'Fore_' harv '_biomes_biom_flux_KKfoodweb_prods.mat'],'bnames',...
    'bios','flux');




