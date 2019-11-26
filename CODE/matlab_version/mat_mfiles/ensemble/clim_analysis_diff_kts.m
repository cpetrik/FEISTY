% Calc AIC of diff k_t values
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/'];


%% Loop over kts
% POEM file info
harv = 'All_fish03';
kays = [0.0405,0.0555,0.0705,0.0855,0.1005,0.1155,0.1305];
skays={'.04','.0555','.07','.0855','.10','.1155','.13'};

fish = NaN*ones(66,length(kays));
lme_Fmcatch = fish;
lme_Pmcatch = fish;
lme_Dmcatch = fish;
lme_AllF = fish;
lme_AllP = fish;
r_all = NaN*ones(5,length(kays));
ss_all = NaN*ones(5,length(kays));
rmse_all = NaN*ones(5,length(kays));
mis_all = NaN*ones(length(kays),45,5);

%%
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

    
    %% LME
    [lme_mcatch,lme_mbio,lme_area] = lme_clim_ensem(sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,lp_my,ld_my);
    
    lme_Fmcatch(:,i) = lme_mcatch(:,1);
    lme_Pmcatch(:,i) = (lme_mcatch(:,2)+lme_mcatch(:,4));
    lme_Dmcatch(:,i) = (lme_mcatch(:,3)+lme_mcatch(:,5));
    
    lme_AllF(:,i) = lme_mbio(:,1)+lme_mbio(:,4);
    lme_AllP(:,i) = lme_mbio(:,2)+lme_mbio(:,5)+lme_mbio(:,7);
    
    %% SAU comparison
    % ADD CALC OF RESIDUALS
    [r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);
    r_all(:,i) = r;
    rmse_all(:,i) = rmse;
    ss_all(:,i) = ss;
    mis_all(i,:,:) = mis;
    
end

%%
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
save([dp 'Climatol_ensemble_diff_kts.mat'],...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','lme_AllF','lme_AllP',...
    'r_all','rmse_all','ss_all','mis_all','kays');
save(['Climatol_ensemble_diff_kts.mat'],...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','lme_AllF','lme_AllP',...
    'r_all','rmse_all','ss_all','mis_all','kays');

