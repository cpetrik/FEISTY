% FEISTY 100 LHS parameter sets
% Only last 12 months of 150 years saved 

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

% SAUP data
load([dp 'SAUP_Stock_top10.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

% param ensemble results
load([dp 'Climatol_ensemble_LHS100.mat']);

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param15_100vals.mat']);

% sfile = sim{M};
% sname = sfile(83:end);

%% SAU comparison
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

l10all = [l10sF(keep);l10sP(keep);l10sD(keep)];
sig = var(l10all);

%% Outputs
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%put residuals of all fn types in one vector
mis_all_fn = mis_all(:,:,2:4);
mis_combo = reshape(mis_all_fn,100,45*3);

mis_fn = mis(:,2:4);
mis_fn = reshape(mis_fn,45*3,1);

%logLike
LL = -135/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);
LL0 = -135/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_fn.^2);
LL_all = [LL0;LL];

%LR test
[h,pval,stat] = lratiotest(LL0,LL,15);

%ALL REJECTED, MUCH WORSE THAN ORIG PARAMS








