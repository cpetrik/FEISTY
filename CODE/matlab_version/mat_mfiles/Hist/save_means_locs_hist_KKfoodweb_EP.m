% Save output of Hist at single locations
% for foodweb diagram
% 145 years, mean of 1951-2000
% Calculates export flux (EP) from det_btm and bottom depth
% Uses NPP - EP for flux to zoop
% Could use zoop prod for flux to zoop

clear all
close all

datap = '/Volumes/FEISTY/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

%load('/Volumes/GFDL/Data/Data_grid_hindcast_NOTflipped.mat');
load('/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Data_grid_hindcast_NOTflipped.mat');
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hist_grid_360x200_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
red = [6,10,12:16];
spots=spots(red)';
shelf = 1:2;
olig = 3:5;
upw = 6:7;
fave = [2;7;4];

dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
sname = 'Historic_';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);
load([dpath sname harv '_locs.mat'])
load([dpath sname 'locs_' harv '_lastyr_sum_mean_biom.mat']);

%% POEM means

mlev = [Flev;Plev;Dlev];
F = squeeze(nansum(all_mean(:,1,:)));
P = squeeze(nansum(all_mean(:,2,:)));
D = squeeze(nansum(all_mean(:,3,:)));
B = squeeze(nansum(all_mean(:,4,:)));
conZ = conZm + conZl;

%% Zoop, det, bent
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'cobalt_det_temp_zoop_npp_means.mat']);
grid = csvread([gpath 'grid_csv.csv']);
ID = grid(:,1);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');

% molN/m2 --> g/m2
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mz_mean_hist = mz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
lz_mean_hist = lz_mean_hist * (106.0/16.0) * 12.01 * 9.0;
% molN/m2/s --> g/m2/d
mzprod_mean_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_mean_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_mean_hist = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_mean_hist = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

z_mean = mz_mean_hist + lz_mean_hist;
z_prod = mzprod_mean_hist + lzprod_mean_hist;

z_mean_grid = z_mean(ID);
z_prod_grid = z_prod(ID);

det_grid = det_mean_hist(ID);

mnpp = npp_mean_hist(ID);

depth_grid = GRD.Z;

%% Estimate EP from Martin curve equation
EP_grid = det_grid ./ ((depth_grid./100).^-0.858);

%%
z_mean_locs = z_mean_grid(ids);
z_prod_locs = z_prod_grid(ids);
det_locs = det_grid(ids);
npp_locs = mnpp(ids);
ep_locs = EP_grid(ids);

%% Table
bios(:,1) = npp_locs(red);
bios(:,2) = z_mean_locs(red);
bios(:,3) = F(red);
bios(:,4) = P(red);
bios(:,5) = B(red);
bios(:,6) = D(red);

flux(:,1) = npp_locs(red) - ep_locs(red); %z_prod_locs(red);
flux(:,2) = conZ(red,1);
flux(:,3) = conF(red,2);
flux(:,4) = det_locs(red);
flux(:,5) = conB(red,3);
flux(:,6) = conP(red,3);
flux(:,7) = conF(red,3);

%% save
rsites = spots;
save([dpath sname harv '_red_locs_biom_flux_KKfoodweb_EP.mat'],'rsites',...
    'bios','flux');




