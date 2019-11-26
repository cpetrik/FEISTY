% FEISTY orig parameter set w/diff kt values
% AIC evaluation 
% Multiply the neg F upwelling LME misfits so they weigh more
% Multiply the P misfits < - log10(2) so they weigh more

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/';

% diff kt's
load([dp 'Climatol_ensemble_diff_kts.mat'],'rmse_all','mis_all','kays',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch');

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);


%% SAU comparison
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

% SAUP data
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');
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

%% Multiply the neg F upwelling LME misfits so they weigh more
up = [3;11;13;27;28;29];
[uboth,uid,kid]=intersect(up,keep);
mis_all_F = mis_all(:,:,2);
negF = mis_all_F(:,kid) < 0;
negF2 = mis_all_F;
negF3 = double(negF);
negF3(negF3==1) = 10;
negF3(negF3==0) = 1;
mis_all_F2 = mis_all_F;
mis_all_F2(:,kid) = mis_all_F(:,kid) .* negF3;

%% Multiply the P misfits < - log10(2) so they weigh more
% used to be log10(5) so they weigh more
mis_all_P = mis_all(:,:,3);
negP = mis_all_P < (-1*log10(2));
negP2 = mis_all_P;
negP3 = double(negP);
negP3(negP3==1) = 3;
negP3(negP3==0) = 1;
mis_all_P2 = mis_all_P .* negP3;

%%
mis_all_D = mis_all(:,:,4);
%put residuals of all fn types in one vector
mis_combo = [mis_all_F2,mis_all_P2,mis_all_D];

%%
mis_fn = mis(:,2:4);
nid = find(mis_fn(:,1) < 0);
unid = intersect(kid,nid);
mis_fn(unid,1) = mis_fn(unid,1) .* 10;

pid = find(mis_fn(:,2) < (-1*log10(2)));
mis_fn(pid,2) = mis_fn(pid,2) .* 3;

mis_fn = reshape(mis_fn,45*3,1);


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
caic_srt2 = nan*ones(length(caic_srt)+1,1);
idc_srt2 = nan*ones(length(caic_srt)+1,1);
idc2 = find(caic_srt<caic);
if (isempty(idc2))
    idc2 = 0;
else
    caic_srt2(idc2) = caic_srt(idc2);
    idc_srt2(idc2) = idc(idc2);
end
caic_srt2(idc2(end)+1) = caic;
caic_srt2((idc2(end)+2):end) = caic_srt((idc2(end)+1):end);

idc_srt2(idc2(end)+1) = 0;
idc_srt2((idc2(end)+2):end) = kays(idc((idc2(end)+1):end));

cdel = caic_srt2 - caic_srt2(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

caicv(:,1) = idc_srt2;
caicv(:,2) = caic_srt2;
caicv(:,3) = cdel;
caicv(:,4) = cw;
cT = array2table(caicv,'VariableNames',{'kt','AIC','delta','weight'});
writetable(cT,[dp 'Orig_params_kts_AIC_Fupneg_mult10_Pneg2_mult3.csv'])


%% vis best maps
for i=1:length(kays)
    kt=kays(i);
    tkfn = num2str(1000+int64(1000*kt));
    cfile = ['Dc_enc70-b200_m4-b175-k',tkfn(2:end),...
        '_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100'];
    fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
    if (i==4)
        load([fpath 'Means_bio_prod_fish_Climatol_All_fish03_',cfile,'.mat']);
    else
        load([fpath 'Climatol_All_fish03.mat']);
    end
    
    %% Maps
    vis_climatol_ensem5(cfile,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,pp);
    
end




