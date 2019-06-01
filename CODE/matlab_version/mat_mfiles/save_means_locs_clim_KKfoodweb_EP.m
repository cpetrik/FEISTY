% Save output of POEM Climatology at single locations
% for foodweb diagram
% 150 years, monthly means saved
% Calculates export flux (EP) from det_btm and bottom depth
% Uses NPP - EP for flux to zoop

clear all
close all

datap = '/Volumes/GFDL/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev','T');
sites = T{:,1};
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
sname = 'Climatol_';
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
cpath = ['/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/'];
load([cpath 'cobalt_zoop_biom_means.mat'],'mz_mean_clim','lz_mean_clim','mzloss_mean_clim','lzloss_mean_clim')
load([cpath 'cobalt_det_biom_means.mat'],'det_mean_clim')

gpath='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([gpath 'clim_npp_Dmeans_Ytot.mat'])
load([gpath 'depth_1deg.mat'])

load(['/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%ESM2.6 in mg C m-2 or mg C m-2 d-1
%from mg C m-2 to g(WW) m-2
% 1e-3 g C in 1 mg C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)

z_mean = (mz_mean_clim + lz_mean_clim) * 1e-3 * 9.0;
z_loss = (mzloss_mean_clim+lzloss_mean_clim) * 1e-3 * 9.0;

z_mean_grid = z_mean(ID);
z_loss_grid = z_loss(ID);

det_grid = det_mean_clim(ID) * 1e-3 * 9.0;

mnpp = npp_mean_clim(ID) * 1e-3 * 9.0;

depth_grid = depth(ID);

%% Estimate EP from Martin curve equation
EP_grid = det_grid ./ ((depth_grid./100).^-0.858);

%%
z_mean_locs = z_mean_grid(ids);
z_loss_locs = z_loss_grid(ids);
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

flux(:,1) = npp_locs(red) - ep_locs(red); %z_loss_locs(red);
flux(:,2) = conZ(red,1);
flux(:,3) = conF(red,2);
flux(:,4) = det_locs(red);
flux(:,5) = conB(red,3);
flux(:,6) = conP(red,3);
flux(:,7) = conF(red,3);

%% save
rsites = sites(red);
save([dpath sname harv '_red_locs_biom_flux_KKfoodweb_EP.mat'],'rsites',...
    'bios','flux');




