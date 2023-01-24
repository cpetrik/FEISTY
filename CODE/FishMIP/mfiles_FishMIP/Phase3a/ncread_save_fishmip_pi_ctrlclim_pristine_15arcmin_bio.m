% FEISTY output at all locations
% PI ctrlclim pristine 1/4 degree

clear all
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

%% Annual means
load([fpath 'Means_PI_ctrlclim_pristine_',cfile,'.mat'],'time',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean')

%% 
% PREFERRED (all units = gWW/m2)
allFB = sf_mean + mf_mean;
allPB = sp_mean + mp_mean + lp_mean;
allDB = sd_mean + md_mean + ld_mean;
allBB = b_mean;


%%
save([fpath 'FishMIP_annual_PI_obsclim_pristine_15arcmin.mat'],'time',...
    'allPB','allDB','allBB','allFB');






