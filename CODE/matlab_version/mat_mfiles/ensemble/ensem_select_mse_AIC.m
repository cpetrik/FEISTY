% FEISTY 100 LHS parameter sets
% Only last 12 months of 150 years saved 

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

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

% SAUP data
load([dp 'SAUP_Stock_top10.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

l10all = [l10sF(keep);l10sP(keep);l10sD(keep)];
%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% Outputs
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%put residuals of all fn types in one vector
mis_all_fn = mis_all(:,:,2:4);
mis_combo = reshape(mis_all_fn,100,45*3);

mis_fn = mis(:,2:4);
mis_fn = reshape(mis_fn,45*3,1);

%% AIC if all models have normally distributed errors
% AIC = n*log(sum(resid?2)/n) + 2*K
aic_all = n * log( sum(mis_combo.^2,2) / n);
aic = n * log( sum(mis_fn.^2) / n);

[aic_srt,idx] = sort(aic_all);
aic_srt2 = [aic;aic_srt];
del = aic_srt2 - aic_srt2(1);
w = exp(-0.5*del) ./ sum( exp(-0.5*del) );

% ALL MUCH MUCH WORSE THAN ORIG PARAMS

aicv(:,1) = [0;idx];
aicv(:,2) = aic_srt2;
aicv(:,3) = del;
aicv(:,4) = w;
T = array2table(aicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(T,[nfile 'LHS_param15_100vals_AIC.csv'])

%% Classic AIC 
% AIC = -2*log(L) + 2*K
% log(L) = (-n/2) * log(2*pi*var) - (1/(2*var)) * sum(resid^2)

%logLike
LL = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);
LL0 = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_fn.^2);
LL_all = [LL0;LL];

caic_all = -2 * LL;
caic = -2 * LL0;

[caic_srt,idc] = sort(caic_all);
caic_srt2 = [caic;caic_srt];
cdel = caic_srt2 - caic_srt2(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

% ALL MUCH MUCH WORSE THAN ORIG PARAMS

caicv(:,1) = [0;idc];
caicv(:,2) = caic_srt2;
caicv(:,3) = cdel;
caicv(:,4) = cw;
cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(cT,[nfile 'LHS_param15_100vals_AIC_classic.csv'])

%% Alt likelihood function
% AIC = -2*log(L) + 2*K
%log(L) = ?n/2 · log(MSE?2)

%logLike
aLL = -n/2 * log(sum(mis_combo.^2,2)/n);
aLL0 = -n/2 * log(sum(mis_fn.^2)/n);

aaic_all = -2*aLL;
aaic = -2*aLL0;

[aaic_srt,ida] = sort(aaic_all);
aaic_srt2 = [aaic;aaic_srt];
adel = aaic_srt2 - aaic_srt2(1);
aw = exp(-0.5*adel) ./ sum( exp(-0.5*adel) );

% GIVES SAME RESULT AS NORM-DISTR EG

aaicv(:,1) = [0;ida];
aaicv(:,2) = aaic_srt2;
aaicv(:,3) = adel;
aaicv(:,4) = aw;
aT = array2table(aaicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
%writetable(aT,[nfile 'LHS_param15_100vals_AIC_alt.csv'])

%% Built in Fn
%logLike LL_all 
k = 15*ones(length(LL),1);
[baic_all,bbic_all] = aicbic(LL,k,n);
[baic,bbic] = aicbic(LL0,15,n);

[baic_srt,idb] = sort(baic_all);
baic_srt2 = [baic;baic_srt];
bdel = baic_srt2 - baic_srt2(1);
bw = exp(-0.5*bdel) ./ sum( exp(-0.5*bdel) );

% GIVES SAME ANSWER AS CLASSIC

baicv(:,1) = [0;idb];
baicv(:,2) = baic_srt2;
baicv(:,3) = bdel;
baicv(:,4) = bw;
bT = array2table(baicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(bT,[nfile 'LHS_param15_100vals_AIC_builtin.csv'])






