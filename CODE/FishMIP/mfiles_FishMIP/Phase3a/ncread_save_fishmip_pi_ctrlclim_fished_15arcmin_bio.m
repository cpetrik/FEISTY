% FEISTY output at all locations
% PI ctrlclim pristine 1/4 degree

clear 
close all

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/QuarterDeg/'];

%% Annual means
load([fpath 'Means_PI_ctrlclim_All_fishobs_' cfile '.mat'],'time',...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean',...
    'mf_my','mp_my','md_my',...
    'lp_my','ld_my')

%% 
% PREFERRED (all units = gWW/m2)
allFB = sf_mean + mf_mean;
allPB = sp_mean + mp_mean + lp_mean;
allDB = sd_mean + md_mean + ld_mean;
allBB = b_mean;

%Total Pelagic Density Catch across Artisanal and Industrial sectors
allFC = mf_my;
allPC = mp_my + lp_my;
allDC = md_my + ld_my;

%%
save([fpath 'FishMIP_annual_PI_obsclim_pristine_15arcmin.mat'],'time',...
    'allFB','allPB','allDB','allBB','allPC','allDC','allFC');






