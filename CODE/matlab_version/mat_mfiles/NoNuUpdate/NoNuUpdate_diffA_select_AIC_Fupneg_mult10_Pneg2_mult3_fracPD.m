% FEISTY No Nu Update parameter tests of A
% FracPD in AIC calc

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%% SAUP data
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

sPD = (Plme_mcatch10) ./ (Plme_mcatch10 + Dlme_mcatch10);

l10all = sPD(keep);
%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% climatol baseline parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Climatology/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

% SAU comparison of baseline params
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);
mis_all_F(1,:,:) = mis(:,2);
mis_all_P(1,:,:) = mis(:,3);
mis_all_D(1,:,:) = mis(:,4);
mis_all_PD(1,:,:) = mis(:,5);

%% SAU comparison of baseline params NoNuUpdate
clear lme_mcatch mis
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
npath = ['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile '/'];
load([npath 'NoNuUpdate_LME_clim_fished_',harv,'_' cfile '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

mis_all_F(2,:,:) = mis(:,2);
mis_all_P(2,:,:) = mis(:,3);
mis_all_D(2,:,:) = mis(:,4);
mis_all_PD(2,:,:) = mis(:,5);

%% SAU comparison of A=0.7 NoNuUpdate
clear lme_mcatch mis
cfile1 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A070_Sm025_nmort1_BE08_noCC_RE00100';
npath1 = ['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile1 '/'];
load([npath1 'NoNuUpdate_LME_clim_fished_',harv,'_' cfile1 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

mis_all_F(3,:,:) = mis(:,2);
mis_all_P(3,:,:) = mis(:,3);
mis_all_D(3,:,:) = mis(:,4);
mis_all_PD(3,:,:) = mis(:,5);

%% SAU comparison of A=0.75 NoNuUpdate
clear lme_mcatch mis
cfile2 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A075_Sm025_nmort1_BE08_noCC_RE00100';
npath2 = ['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile2 '/'];
load([npath2 'NoNuUpdate_LME_clim_fished_',harv,'_' cfile2 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

mis_all_F(4,:,:) = mis(:,2);
mis_all_P(4,:,:) = mis(:,3);
mis_all_D(4,:,:) = mis(:,4);
mis_all_PD(4,:,:) = mis(:,5);

%% SAU comparison of A=0.8 params NoNuUpdate
clear lme_mcatch mis
cfile3 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A080_Sm025_nmort1_BE08_noCC_RE00100';
npath3 = ['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile3 '/'];
load([npath3 'NoNuUpdate_LME_clim_fished_',harv,'_' cfile3 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

mis_all_F(5,:,:) = mis(:,2);
mis_all_P(5,:,:) = mis(:,3);
mis_all_D(5,:,:) = mis(:,4);
mis_all_PD(5,:,:) = mis(:,5);

%% put residuals of all fn types in one vector
mis_combo = mis_all_PD;

%id = {'Orig','A50','A70','A75','A80'};
id = [0,50,70,75,80];

%% Classic AIC 
% AIC = -2*log(L) + 2*K
% log(L) = (-n/2) * log(2*pi*var) - (1/(2*var)) * sum(resid^2)

%logLike
LL = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);

caic_all = -2 * LL;

[caic_srt,idc] = sort(caic_all);
id_srt = id(idc);

%%
cdel = caic_srt - caic_srt(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

%caicv(:,1) = idc;
caicv(:,1) = id(idc);
caicv(:,2) = caic_srt;
caicv(:,3) = cdel;
caicv(:,4) = cw;
%cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
cT = array2table(caicv,'VariableNames',{'A','AIC','delta','weight'});
writetable(cT,'/Volumes/FEISTY/NC/Clim_comp_tests/NoNuUpdate_diffA_AIC_FracPD.csv')

