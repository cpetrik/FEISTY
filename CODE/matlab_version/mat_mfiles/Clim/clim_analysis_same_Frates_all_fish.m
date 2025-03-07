% FEISTY F rate search for MMSSY
% Same F rate applied to all fish
% Climatology

clear all
close all

load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00001/'];

Fish = [0.0:0.1:2 2.2:0.2:5];

%%
SF = NaN*ones(45124,length(Fish));
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
lme_Fmcatch = NaN*ones(66,length(Fish));
lme_Pmcatch = NaN*ones(66,length(Fish));
lme_Dmcatch = NaN*ones(66,length(Fish));
r_all = NaN*ones(5,length(Fish));
ss_all = NaN*ones(5,length(Fish));
rmse_all = NaN*ones(5,length(Fish));
mis_all = NaN*ones(length(Fish),45,5);

%%

for M=1:length(Fish)
    frate = Fish(M);
    tfish = num2str(100+int64(10*frate));
    if (frate > 0)
        sfile = [dp 'Climatol_All_fish' tfish(2:end)];
    else
        sfile = [dp 'Climatol_pristine'];
    end
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
    
    %% LME
    [lme_mcatch,lme_mbio,lme_area] = lme_clim_ensem(sf_mean,sp_mean,sd_mean,...
        mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,lp_my,ld_my);
    
    %% Save
    save([sfile '.mat'],...
        'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
        'lp_mean','ld_mean','time','lyr',...
        'mf_my','mp_my','md_my','lp_my','ld_my',...
        'lme_mcatch','lme_mbio','lme_area','-append');
    
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
    
    %% SAU comparison
    % ADD CALC OF RESIDUALS
    [r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);
    r_all(:,M) = r;
    rmse_all(:,M) = rmse;
    ss_all(:,M) = ss;
    mis_all(M,:,:) = mis;
    
end

%%
save([dp 'Climatol_same_Frates_all_fish.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch',...
    'r_all','rmse_all','ss_all','mis_all','Fish')

