%Compare BGC offline and online

clear
close all

%% Offline
offp = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_only/';

%update with results from my run next week
load([offp '20040101.ocean_feisty_forcing_means_10yr_remy.mat'])
load([offp '20040101.ocean_cobalt_10yr_BATS_remy.mat'])
load([offp '20040101.ocean_cobalt_tracers_daily_means_10yr_remy.mat'])
load([offp '20040101.ocean_grid_12mo_BATS.mat'])

%% Online
onp = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/cobalt_feisty/test_BATS_10yr/';

load([offp '20040101.ocean_feisty_forcing_means_10yr.mat'])
load([offp '20040101.ocean_cobalt_tracers_month_z_means.mat'])

%%
