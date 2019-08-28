% FEISTY full parameter sets of 5 most sensitive params
% Only high, mid, and low values 
% Remove all F outliers from analysis 
% Need to weight F points so equal to P points b/c numbers differ?

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

% param ensemble results
load([dp 'Climatol_ensemble_param5_mid5.mat'],'rmse_all','mis_all','sim',...
   'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch');
% load(['Climatol_ensemble_param5_mid5.mat'],'rmse_all','mis_all','sim',...
%     'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch');

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param5_mid5.mat']);

% sfile = sim{M};
% sname = sfile(83:end);

%% SAU comparison
[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);

% SAUP data
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;
%Drop LMEs with overestimate of F catch (log10SAU < -2)
over = [1;2;6;7;12;18;53;54;65]; 
[oboth,kid]=setdiff(keep,over);

l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

l10all = [l10sF(oboth);l10sP(keep);l10sD(keep)];
%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% Outputs
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%put residuals of all fn types in one vector
mis_all_F = mis_all(:,kid,2);   %36 F obs
mis_all_PD = mis_all(:,:,3:4);  %45 P&D obs
mis_PD = reshape(mis_all_PD,length(fx_all),45*2);
%mult F obs by difference in # of obs, so have the same "weight"
%then mult by extra 2 because F still too low
mis_all_F2 = mis_all_F * (45/36)*2;
mis_combo = [mis_all_F2,mis_PD];

mis_fn = mis(:,2:4);
mis_fn = reshape(mis_fn,45*3,1);
mis_fn = mis_fn(~isnan(mis_fn));

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
idc2 = find(caic_srt<caic);

caic_srt2 = nan*ones(length(caic_srt)+1,1);
caic_srt2(idc2) = caic_srt(idc2);
caic_srt2(idc2(end)+1) = caic;
caic_srt2((idc2(end)+2):end) = caic_srt((idc2(end)+1):end);

idc_srt2 = nan*ones(length(caic_srt)+1,1);
idc_srt2(idc2) = idc(idc2);
idc_srt2(idc2(end)+1) = 0;
idc_srt2((idc2(end)+2):end) = idc((idc2(end)+1):end);

cdel = caic_srt2 - caic_srt2(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

caicv(:,1) = idc_srt2;
caicv(:,2) = caic_srt2;
caicv(:,3) = cdel;
caicv(:,4) = cw;
cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(cT,[dp 'LHS_param5_mid5_AIC_Foutlier_weighted2.csv'])

%% Built in Fn
%logLike LL_all 
k = 5*ones(length(LL),1);
[baic_all,bbic_all] = aicbic(LL,k,n);
[baic,bbic] = aicbic(LL0,5,n);

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
writetable(bT,[dp 'LHS_param5_mid5_AIC_builtin_Foutlier_weighted2.csv'])


%% AICs <= AIC(orig) + 2
besti = find(baicv(:,2) <= baic+2);
pid = baicv(besti,1);
pid = pid(pid>0);
pset(:,1) = pid;
pset(:,2:6) = fx_all(pid,:);
pset(:,7) = baic_all(pid);

pT = array2table(pset,'VariableNames',{'ParamSet','Lambda','bMet','bEnc',...
    'aMet','aEnc','AIC'});
writetable(pT,[dp 'LHS_param5_mid5_bestAIC_params_Foutlier_weighted2.csv'])

id1 = pid;

% NOTHING NEARLY AS GOOD AS ORIGINAL PARAMETER SET

%% vis best maps
nsim = length(id1);
for j=1:nsim
    M=id1(j);
    sfile = sim{M};
    sname = sfile(140:end);
    load([sfile '_means.mat']);
    
    %% Maps
    vis_climatol_ensem5(sname,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,pp);
    
end

% P & D ARE BAD IF A_MET = 6

%% vis best SAUP comp
for j=1:length(id1)
    M=id1(j);
    sfile = sim{M};
    sname = sfile(140:end);
    
    %% Comp
    vis_clim_lme_saup_corr_stock(lme_Fmcatch(:,M),lme_Pmcatch(:,M),...
        lme_Dmcatch(:,M),pp,sname);
    
end

