% FEISTY full parameter sets of 5 most sensitive params
% Only high, mid, and low values
% Compare:  multFup_neg5_multPneg3
%           multFup_neg2_multPneg3
%           multFup_neg2_multPneg2

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

%%
load([dp 'LHS_param5_mid5_bestAIC_params_multFup_neg_multPneg.mat']);
params1 = params;
pid1 = pid;
clear params pid

load([dp 'LHS_param5_mid5_bestAIC_params_multFup_neg2_multPneg3.mat']);
params3 = params;
pid3 = pid;
clear params pid

load([dp 'LHS_param5_mid5_bestAIC_params_multFup_neg2_multPneg2.mat']);
params2 = params;
pid2 = pid;
clear params pid

%%
nd2 = intersect(pid1,pid2);
nd3 = intersect(pid1,pid3);

dd2 = setdiff(pid2,pid1);

