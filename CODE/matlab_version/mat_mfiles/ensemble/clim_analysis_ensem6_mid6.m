% FEISTY full parameter sets of 5 most sensitive params
% Only last 12 months of 150 years saved 
% Only high, mid, and low values

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
dp = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

%%
nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param6_mid6_temp.mat']);

sim = cell(length(fx),1);
for j = 1:length(fx)
    %! Change individual parameters
    pset = fx(j,:);
    %pset=round(pset,3);
    set_params6(pset)
    
    %! Make core parameters/constants (global)
    const_params6()

    %! Create a directory for output
    fname = sub_fname_ensemble6();
    sim{j} = fname;
end

%%
SF = NaN*ones(45124,length(fx));
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
lme_Fmcatch = NaN*ones(66,length(fx));
lme_Pmcatch = NaN*ones(66,length(fx));
lme_Dmcatch = NaN*ones(66,length(fx));
lme_AllF = NaN*ones(66,length(fx));
lme_AllP = NaN*ones(66,length(fx));
r_all = NaN*ones(5,length(fx));
ss_all = NaN*ones(5,length(fx));
rmse_all = NaN*ones(5,length(fx));
mis_all = NaN*ones(length(fx),45,5);

%%
for M=1:length(fx)
    sfile = sim{M};
    sname = sfile(131:end);
    load([sfile '.mat']);
    
    %% Last year means
    [id,nt] = size(Spinup_Bent.bio);
    time=1:nt;
    lyr=time((end-12+1):end);
    sp_mean=mean(Spinup_Sml_p.bio(:,lyr),2);
    sf_mean=mean(Spinup_Sml_f.bio(:,lyr),2);
    sd_mean=mean(Spinup_Sml_d.bio(:,lyr),2);
    mp_mean=mean(Spinup_Med_p.bio(:,lyr),2);
    mf_mean=mean(Spinup_Med_f.bio(:,lyr),2);
    md_mean=mean(Spinup_Med_d.bio(:,lyr),2);
    lp_mean=mean(Spinup_Lrg_p.bio(:,lyr),2);
    ld_mean=mean(Spinup_Lrg_d.bio(:,lyr),2);
    b_mean=mean(Spinup_Bent.bio(:,lyr),2);
    
    mf_my=mean(Spinup_Med_f.yield(:,lyr),2);
    mp_my=mean(Spinup_Med_p.yield(:,lyr),2);
    md_my=mean(Spinup_Med_d.yield(:,lyr),2);
    lp_my=mean(Spinup_Lrg_p.yield(:,lyr),2);
    ld_my=mean(Spinup_Lrg_d.yield(:,lyr),2);
    
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
    
    %% Maps
%     vis_climatol_ensem(sname,sf_mean,sp_mean,sd_mean,...
%     mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean);
    
    %% LME
    [lme_mcatch,lme_mbio,lme_area] = lme_clim_ensem(sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,lp_my,ld_my);
    
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
    
    %% Save
    save([sfile '_means.mat'],...
        'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
        'lp_mean','ld_mean','time','lyr',...
        'mf_my','mp_my','md_my','lp_my','ld_my',...
        'lme_mcatch','lme_mbio','lme_area');
    
end

%%
save([dp 'Climatol_ensemble_param6_mid6_temp.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','lme_AllF','lme_AllP',...
    'r_all','rmse_all','ss_all','mis_all','sim','fx')

