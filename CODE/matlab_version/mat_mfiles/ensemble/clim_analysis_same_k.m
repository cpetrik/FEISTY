% FEISTY climatology/baseline parameter set 
% kE = kC = kM (kt) for all

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/'...
'Dc_enc70-b200_m4-b175_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100_kappaA50/'];

%%
kays = [63,76,86,96];

SF = NaN*ones(45124,length(kays));
SP = SF;
SD = SF;
MF = SF;
MP = SF;
MD = SF;
LP = SF;
LD = SF;
BI = SF;
MFc = SF;
MPc = SF;
MDc = SF;
LPc = SF;
LDc = SF;
lme_Fmcatch = NaN*ones(66,length(kays));
lme_Pmcatch = NaN*ones(66,length(kays));
lme_Dmcatch = NaN*ones(66,length(kays));
lme_AllF = NaN*ones(66,length(kays));
lme_AllP = NaN*ones(66,length(kays));
r_all = NaN*ones(5,length(kays));
ss_all = NaN*ones(5,length(kays));
rmse_all = NaN*ones(5,length(kays));
mis_all = NaN*ones(length(kays),45,5);

%%
sim = cell(length(kays),1);
for M=1:length(kays)
    sname = ['Dc_enc70-b200-k0',num2str(kays(M)),'_m4-b175-k0',...
        num2str(kays(M)),'_c20-b250-k0',num2str(kays(M)),...
        '_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100_kappaA50'];
    sfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' sname];
    
    load([sfile '/Means_Climatol_All_fish03_' sname '.mat']);
    
    sim{M} = sname;
    
    %% Maps
    %     vis_climatol_ensem(sname,sf_mean,sp_mean,sd_mean,...
    %     mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean);
    
    %% LME
    [lme_mcatch,lme_mbio,lme_area] = lme_clim_ensem(sf_mean,sp_mean,sd_mean,...
        mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,lp_my,ld_my);
    
    %% Save
    save([sfile '/LME_means_Climatol_All_fish03_' sname '.mat'],...
        'lme_mcatch','lme_mbio','lme_area');
    
    %% All sims
    SF(:,M) = sf_mean;
    SP(:,M) = sp_mean;
    SD(:,M) = sd_mean;
    MF(:,M) = mf_mean;
    MP(:,M) = mp_mean;
    MD(:,M) = md_mean;
    LP(:,M) = lp_mean;
    LD(:,M) = ld_mean;
    BI(:,M) = b_mean;
    
    MFc(:,M) = mf_my;
    MPc(:,M) = mp_my;
    MDc(:,M) = md_my;
    LPc(:,M) = lp_my;
    LDc(:,M) = ld_my;
    
    lme_Fmcatch(:,M) = lme_mcatch(:,1);
    lme_Pmcatch(:,M) = (lme_mcatch(:,2)+lme_mcatch(:,4));
    lme_Dmcatch(:,M) = (lme_mcatch(:,3)+lme_mcatch(:,5));
    
    lme_AllF(:,M) = lme_mbio(:,1)+lme_mbio(:,4);
    lme_AllP(:,M) = lme_mbio(:,2)+lme_mbio(:,5)+lme_mbio(:,7);
    
    %% SAU comparison
    % ADD CALC OF RESIDUALS
    [r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);
    r_all(:,M) = r;
    rmse_all(:,M) = rmse;
    ss_all(:,M) = ss;
    mis_all(M,:,:) = mis;
    
end

%%
save([dp 'Climatol_baseline_same_ks.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','lme_AllF','lme_AllP',...
    'r_all','rmse_all','ss_all','mis_all','sim','kays')

