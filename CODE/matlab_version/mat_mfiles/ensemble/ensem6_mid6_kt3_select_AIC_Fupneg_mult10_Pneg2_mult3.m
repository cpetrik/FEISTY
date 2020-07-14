% FEISTY full parameter sets of 6 most sensitive params
% Only high, mid, and low values

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp  = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
dp2 = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

%% param ensemble results with kt=0.0885
load([dp 'Climatol_ensemble_param5_mid5.mat'],'rmse_all','mis_all','sim',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','r_all');
r_all1 = r_all;
rmse_all1 = rmse_all;
mis_all1 = mis_all;
sim1 = sim;
lme_Fmcatch1 = lme_Fmcatch;
lme_Pmcatch1 = lme_Pmcatch;
lme_Dmcatch1 = lme_Dmcatch;
clear rmse_all mis_all sim lme_Fmcatch lme_Pmcatch lme_Dmcatch r_all

% param ensemble results with kt=0.0655 and 0.0755
load([dp2 'Climatol_ensemble_param6_mid6_temp2.mat'],'rmse_all','mis_all','sim',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','r_all');
r_all2 = r_all;
rmse_all2 = rmse_all;
mis_all2 = mis_all;
sim2 = sim;
lme_Fmcatch2 = lme_Fmcatch;
lme_Pmcatch2 = lme_Pmcatch;
lme_Dmcatch2 = lme_Dmcatch;
clear rmse_all mis_all sim lme_Fmcatch lme_Pmcatch lme_Dmcatch r_all

% param ensemble results with kt=0.0955 
load([dp2 'Climatol_ensemble_param6_mid6_temp3.mat'],'rmse_all','mis_all','sim',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','r_all');
r_all3 = r_all;
rmse_all3 = rmse_all;
mis_all3 = mis_all;
sim3 = sim;
lme_Fmcatch3 = lme_Fmcatch;
lme_Pmcatch3 = lme_Pmcatch;
lme_Dmcatch3 = lme_Dmcatch;
clear rmse_all mis_all sim lme_Fmcatch lme_Pmcatch lme_Dmcatch r_all

%%
r_all = [r_all1, r_all2, r_all3];
rmse_all = [rmse_all1, rmse_all2, rmse_all3];
mis_all = [mis_all1; mis_all2; mis_all3];
sim = [sim1; sim2; sim3];
lme_Fmcatch = [lme_Fmcatch1, lme_Fmcatch2, lme_Fmcatch3];
lme_Pmcatch = [lme_Pmcatch1, lme_Pmcatch2, lme_Pmcatch3];
lme_Dmcatch = [lme_Dmcatch1, lme_Dmcatch2, lme_Dmcatch3];

% climatol parameter set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Climatology/'];
load([dpath 'LME_clim_fished_',harv,'_' cfile '.mat']);

%%
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param5_mid5.mat'],'fx_all');
load([nfile 'LHS_param6_mid6_temp2.mat'],'fx','ptext');
load([nfile 'LHS_param6_mid6_temp3.mat'],'fz3');

ptext{6}='kt';

fx_all(:,6) = 0.0855;
fx_all = [fx_all; fx; fz3];

% sfile = sim{M};
% sname = sfile(83:end);

save([dp2 'Climatol_ensemble_param6_mid6_temp6789.mat'],'r_all','rmse_all',...
    'mis_all','sim','lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','fx_all');

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
writetable(cT,[dp 'LHS_param6_mid6_kt2_AIC_Fupneg_mult10_Pneg2_mult3.csv'])

%% Built in Fn - DON'T HAVE THIS TOOLBOX ANYMORE
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
% writetable(bT,[dp 'LHS_param6_mid6_AIC_builtin_Fupneg_mult10_Pneg2_mult3.csv'])


%% AICs <= AIC(orig) + 2
besti = find(caicv(:,2) <= caic+2);
pid = caicv(besti,1);
pid = pid(pid>0);
pset(:,1) = pid;
pset(:,2:7) = fx_all(pid,:);
pset(:,8) = caic_all(pid);
pset(:,9:12) = rmse_all(1:4,pid)';
pset(:,13:16) = r_all(1:4,pid)';

pT = array2table(pset,'VariableNames',{'ParamSet','Lambda','bMet','bEnc',...
    'aMet','aEnc','kMet','AIC','rmseAll','rmseF','rmseP','rmseD','rAll',...
    'rF','rP','rD'});
writetable(pT,[dp 'LHS_param6_mid6_kt2_bestAIC_params_Fupneg_mult10_Pneg2_mult3.csv'])

id1 = pid;

params = fx_all(pid,:);
save([dp 'LHS_param6_mid6_kt2_bestAIC_params_Fupneg_mult10_Pneg2_mult3.mat'],...
    'params','ptext','pid')

% 51 Sets (added 17 sets with kt=0.0955)
test=find(pset(:,7)==0.0955);

%% vis best maps
% 8 sets don't have great looking F
bad = id1([27,33,35,36,40,41,48,50]);
% bad = [429;248;409;578;427;411;471;407]; same as before
id1 = setdiff(pid,bad);
for j=1:length(id1)
    M=id1(j);
    
    if (M >=730)
        sname = sim{M};
        sfile = [dp2 sname];
    elseif (M <= 243)
        sfile = sim{M};
        sname = sfile(140:end); 
    else
        sfile = sim{M};
        sname = sfile(131:end);
    end
    load([sfile '_means.mat']);
    
    %% Maps
    vis_climatol_ensem5(sname,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,pp);
    
end

% 8 sets don't have great looking F - so 43 good sets
% 7/8 of these sets have the lowest D rmse values
% [429;248;409;578;427;411;471;407]
red_params = fx_all(id1,:);
save([dp 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext','id1')

pT = array2table(red_params,'VariableNames',ptext);
writetable(pT,[dp 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.csv'])

[blah,rid] = intersect(pset(:,1),id1);
rpset = pset(rid,1:8);
rT = array2table(rpset,'VariableNames',{'ParamSet','Lambda','bMet','bEnc',...
    'aMet','aEnc','kMet','AIC'});
writetable(rT,[dp 'LHS_param6_mid6_kt2_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced_AIC.csv'])

%% vis best SAUP comp
for j=1:length(id1)
    M=id1(j);
    sfile = sim{M};
    sname = sfile(140:end);
    
    %% Comp
    vis_clim_lme_saup_corr_stock(lme_Fmcatch(:,M),lme_Pmcatch(:,M),...
        lme_Dmcatch(:,M),pp,sname);
    
end


