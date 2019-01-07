% POEM output at all locations
% Only last 12 months of 150 years saved 
% +/-10% changes to parameters

clear all
close all

ptext = {'base','h45','h55','gam45','gam55','amet036','amet044','lam063',...
    'lam077','bc189','bc231','be189','be231','bm189','bm231','BE0675','BE0825',...
    'RE0009','RE0011','fish027','fish033','kc0544','kc0664','ke0544','ke0664',...
    'kt0544','kt0664','kap45','kap55','unat09','unat11','A063','A077',...
    'Sm0225','Sm0275','J063','J077','D063','D077'};

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
    
SF = NaN*ones(45124,length(ptext));
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

for M=1:length(ptext)
    pt = ptext{M};
    
    sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/param_sens/Climatol_All_fish03_',pt,'.mat'];
    load(sfile);
    
    % Last year
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
    
    %
    save(sfile,...
        'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
        'lp_mean','ld_mean','time','lyr',...
        'mf_my','mp_my','md_my','lp_my','ld_my',...
        '-append');
    
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
    
end

%%
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
save([nfile 'Climatol_All_fish03_means_param_sens_10p.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc','ptext')

%% Outputs
% Relative biomass of each group
% Total biomass
% Weighted latitudinal center of biomass





