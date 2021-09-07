% FEISTY F rate search for MMSSY
% Same F rate applied to all fish
% Loop through F and RE pairs
% Climatology

clear all
close all

load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00001/'];

Fi = logspace(-1,0.7,11);
Fish = Fi([1,3:end]);
reff = logspace(-3,-1,10);
reff = reff(2:end);

%%
SF = NaN*ones(length(Fish),length(reff),45124);
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
lme_Fmcatch = NaN*ones(length(Fish),length(reff),66);
lme_Pmcatch = lme_Fmcatch;
lme_Dmcatch = lme_Fmcatch;
r_all = NaN*ones(length(Fish),length(reff),5);
rmse_all = NaN*ones(length(Fish),length(reff),5);

%%
for k=1:length(reff)
    for M=1:length(Fish)
        frate = Fish(M);
        rfrac = reff(k);
        
        tfish = num2str(100+int64(10*frate));
        tre = num2str(100000+int64(round(10000*rfrac)));
        
        dp = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
        'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_',...
        'noCC_RE',tre(2:end),'/'];

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
        SF(M,k,:) = sf_mean;
        SP(M,k,:) = sp_mean;
        SD(M,k,:) = sd_mean;
        MF(M,k,:) = mf_mean;
        MP(M,k,:) = mp_mean;
        MD(M,k,:) = md_mean;
        LP(M,k,:) = lp_mean;
        LD(M,k,:) = ld_mean;
        BI(M,k,:) = b_mean;
        
        MFc(M,k,:) = mf_my;
        MPc(M,k,:) = mp_my;
        MDc(M,k,:) = md_my;
        LPc(M,k,:) = lp_my;
        LDc(M,k,:) = ld_my;
        
        lme_Fmcatch(M,k,:) = lme_mcatch(:,1);
        lme_Pmcatch(M,k,:) = (lme_mcatch(:,2)+lme_mcatch(:,4));
        lme_Dmcatch(M,k,:) = (lme_mcatch(:,3)+lme_mcatch(:,5));
        
        %% SAU comparison
        % ADD CALC OF RESIDUALS
        [r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch);
        r_all(M,k,:) = r;
        rmse_all(M,k,:) = rmse;
        
    end
end

%%
save([dp 'Climatol_same_Frates_all_fish_RE_sens.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch',...
    'r_all','rmse_all','Fish','reff')

