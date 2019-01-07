% Catch and P:D ratio by LME 
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results & SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cdir 'temp_100_1deg_ESM26_5yr_clim_191_195.mat'])
load([cdir 'btm_temp_1deg_ESM26_5yr_clim_191_195.mat'])
load([cpath 'LME_clim_temp_zoop_det.mat']);

ptemp_mean_clim=squeeze(nanmean(temp_100,1));
btemp_mean_clim=squeeze(nanmean(btm_temp,1));
tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

%% POEM
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/'];
dpath = [dp cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_area','rPD_biom','rPD_catch','rPD_catch_mtkm2');

lme_area_km2 = lme_area * 1e-6;

%% DvD 
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
%FracLP(L)
% all LMEs
LME_all = [1:66]';
FracLP_all = FracLP;
FracLP_all(64,1) = NaN;
FracLP_all(65,1) = NaN;
FracLP_all(66,1) = NaN;

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

%% POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;

% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);

%% Add small num to 0s / subtract small num from 1s
D0=find(FracLP_all==0);
D1=find(FracLP_all==1);
FracLP_all(D0) = FracLP_all(D0)+eps;
FracLP_all(D1) = FracLP_all(D1)-1.1e-6;

%% Data for comparison stats
Fstat = table(LME_all,l10s,l10sF,l10sP,l10sD,...
    l10p,l10pF,l10pP,l10pD,sFracPD,pFracPD,FracLP_all,...
    'VariableNames',{'LME','SAll','SF','SP','SD','PAll','PF','PP','PD',...
    'SPD','PPD','DPD'});
writetable(Fstat,[dpath 'LME_DvD_SAU_vals_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)

Fstat2 = table(keep,l10s(keep),l10sF(keep),l10sP(keep),l10sD(keep),...
    l10p(keep),l10pF(keep),l10pP(keep),l10pD(keep),...
    sFracPD(keep),pFracPD(keep),FracLP_all(keep),...
    'VariableNames',{'LME','SAll','SF','SP','SD','PAll','PF','PP','PD',...
    'SPD','PPD','DPD'});
writetable(Fstat2,[dpath 'LME_DvD_SAU_notLELC_vals_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)

save([dpath 'LME_DvD_SAU_vals_' cfile '.mat'],'Fstat','Fstat2')

