% Save output of POEM Climatology at single locations
% for foodweb diagram
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

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

% Colors
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([3,1,5],:);
cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

mass = [M_s;M_m;M_l];
mass = repmat(mass,1,length(spots));
L = [L_s;L_m;L_l];

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

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

z_mean_locs = z_mean_grid(ids);
z_loss_locs = z_loss_grid(ids);
det_locs = det_grid(ids);
npp_locs = mnpp(ids);

%% Table
bios(:,1) = npp_locs(red);
bios(:,2) = z_mean_locs(red);
bios(:,3) = F(red);
bios(:,4) = P(red);
bios(:,5) = B(red);
bios(:,6) = D(red);

flux(:,1) = z_loss_locs(red);
flux(:,2) = conZ(red,1);
flux(:,3) = conF(red,2);
flux(:,4) = det_locs(red);
flux(:,5) = conB(red,3);
flux(:,6) = conP(red,3);
flux(:,7) = conF(red,3);

%% save
rsites = sites(red);
save([dpath sname harv '_red_locs_biom_flux_KKfoodweb.mat'],'rsites',...
    'bios','flux');




