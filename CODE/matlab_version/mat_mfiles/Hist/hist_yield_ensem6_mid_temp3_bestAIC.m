% FEISTY Historic runs of best parameter sets
% varying 6 most sensitive params (added kt)

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';

pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Hist_param6_mid_best/'];
if (~isfolder(pp))
    mkdir(pp)
end

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
params = red_params;

epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Hlme_mcatch')

%% Hist grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','AREA_OCN');
grid = csvread([gpath 'grid_csv.csv']);
load([gpath 'lme_mask_esm2m.mat']);
load([cpath 'LME_hist9095_temp_zoop_det.mat'],'lme_ptemp','lme_area');

hlme_ptemp = lme_ptemp;
clear lme_ptemp

hlme_area_km2 = lme_area * 1e-6;
clear lme_area

hlme = lme_mask_esm2m';
hID = grid(:,1);

%% Catch to compare to SAU
hlme_mcatch = squeeze(nansum(Hlme_mcatch,2)) * 1e-6;
hlme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
hlme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
hlme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;

% MT/km2
hlme_mcatch_MTkm2 = hlme_mcatch ./ repmat(hlme_area_km2,1,43);
hlme_Fmcatch_MTkm2 = hlme_Fmcatch ./ hlme_area_km2;
hlme_Pmcatch_MTkm2 = hlme_Pmcatch ./ hlme_area_km2;
hlme_Dmcatch_MTkm2 = hlme_Dmcatch ./ hlme_area_km2;

% totals
hlme_tcatch_MTkm2 = sum(hlme_mcatch_MTkm2);
hlme_Ftcatch_MTkm2 = sum(hlme_Fmcatch_MTkm2);
hlme_Ptcatch_MTkm2 = sum(hlme_Pmcatch_MTkm2);
hlme_Dtcatch_MTkm2 = sum(hlme_Dmcatch_MTkm2);

%% SAUP: All, F, P, D, Frac P
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

