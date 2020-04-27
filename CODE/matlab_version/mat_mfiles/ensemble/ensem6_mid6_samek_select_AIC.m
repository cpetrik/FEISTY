% FEISTY full parameter sets of 5 most sensitive params + k
% Only high, mid, and low values
% vary all k's, set as same value
% AIC with equal weighting

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/',...
    'param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/',...
    'bestAIC_params/'];

% param ensemble results
load([dp 'Climatol_ensemble_param6_mid6_samek.mat'],'rmse_all','mis_all',...
    'r_all','fsim','sim','lme_Fmcatch','lme_Pmcatch','lme_Dmcatch');

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param6_mid6_samek.mat'],'fx');

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

%put residuals of all fn types in one vector
mis_all_fn = mis_all(:,:,2:4);
mis_combo = reshape(mis_all_fn,length(fx),45*3);

mis_fn = mis(:,2:4);
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
writetable(cT,[nfile 'LHS_param6_mid6_samek_AIC_classic.csv'])

%% Built in Fn - in Economterics toolbox, no longer have
%logLike LL_all 
% k = 5*ones(length(LL),1);
% [baic_all,bbic_all] = aicbic(LL,k,n);
% [baic,bbic] = aicbic(LL0,5,n);
% 
% [baic_srt,idb] = sort(baic_all);
% baic_srt2 = [baic;baic_srt];
% bdel = baic_srt2 - baic_srt2(1);
% bw = exp(-0.5*bdel) ./ sum( exp(-0.5*bdel) );
% 
% % GIVES SAME ANSWER AS CLASSIC
% 
% baicv(:,1) = [0;idb];
% baicv(:,2) = baic_srt2;
% baicv(:,3) = bdel;
% baicv(:,4) = bw;
% bT = array2table(baicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
% writetable(bT,[nfile 'LHS_param6_mid6_samek_AIC_builtin.csv'])


%% AICs <= AIC(orig) + 2
besti = find(caicv(:,2) <= caic+2);
pid = caicv(besti,1);
pid = pid(pid>0);
pset(:,1) = pid;
pset(:,2:7) = fx(pid,:);
pset(:,8) = caic_all(pid);
pset(:,9:12) = rmse_all(1:4,pid)';
pset(:,13:16) = r_all(1:4,pid)';

pT = array2table(pset,'VariableNames',{'ParamSet','Lambda','bMet','bEnc',...
    'aMet','aEnc','kMet','AIC','rmseAll','rmseF','rmseP','rmseD','rAll',...
    'rF','rP','rD'});
writetable(pT,[nfile 'LHS_param6_mid6_samek_bestAIC_params.csv'])

id1 = pid; %78 sets!

%% Plot RMSE in 3D space
figure(1)
subplot(2,2,1)
scatter3(rmse_all(3,:),rmse_all(4,:),rmse_all(2,:)); hold on;
scatter3(rmse_all(3,id1),rmse_all(4,id1),rmse_all(2,id1),'r','filled'); hold on;
scatter3(rmse(3),rmse(4),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
zlabel('F RMSE')
% xlim([0.25 8])
% ylim([0.25 1])
% zlim([0 6])
subplot(2,2,2)
scatter(rmse_all(3,:),rmse_all(4,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(4,id1),'r','filled'); hold on;
scatter(rmse(3),rmse(4),'k','filled');
xlabel('P RMSE')
ylabel('D RMSE')
% xlim([0.25 8])
% ylim([0.25 1])
subplot(2,2,3)
scatter(rmse_all(3,:),rmse_all(2,:)); hold on;
scatter(rmse_all(3,id1),rmse_all(2,id1),'r','filled'); hold on;
scatter(rmse(3),rmse(2),'k','filled');
xlabel('P RMSE')
ylabel('F RMSE')
% xlim([0.25 8])
% ylim([0 6])
subplot(2,2,4)
scatter(rmse_all(4,:),rmse_all(2,:)); hold on;
scatter(rmse_all(4,id1),rmse_all(2,id1),'r','filled'); hold on;
scatter(rmse(4),rmse(2),'k','filled');
xlabel('D RMSE')
ylabel('F RMSE')
% xlim([0.25 1])
% ylim([0 6])
print('-dpng',[pp 'RMSE_SAUP_type_best_AIC_param6_mid6_samek.png'])

%% vis best maps
for j=1:length(id1)
    M=id1(j);
    sfile = fsim{M};
    sname = sfile(140:end);
    load([sfile '_means.mat']);
    
    %% Maps
    vis_climatol_ensem5(sname,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,pp);
    
end

% Over half have no F in tropics/subtropics

%% vis best SAUP comp
for j=1:length(id1)
    M=id1(j);
    sfile = sim{M};
    sname = sfile(140:end);
    
    %% Comp
    vis_clim_lme_saup_corr_stock(lme_Fmcatch(:,M),lme_Pmcatch(:,M),...
        lme_Dmcatch(:,M),pp,sname);
    
end







