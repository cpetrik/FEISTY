% POEM output at all locations
% Only last 12 months of 100 years saved to compare to SAUP for tuning

clear all
close all

%kays = [0.0405,0.0555,0.0705,0.1005,0.1155,0.1305];
kays = [0.025,0.05,0.1,0.125];

for M=1:length(kays)
    kt = 0.0855;
    bent_eff = kays(M);
    tbe = num2str(100+int64(100*bent_eff));
    tkfn = num2str(1000+int64(1000*kt));
    
    cfile = ['Dc_enc70-b200_m4-b175-k',tkfn(2:end),...
        '_c20-b250_D075_J100_A050_Sm025_nmort1_BE',tbe(2:end),...
        '_noCC_RE00100'];
    sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/Climatol_All_fish03.mat'];
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
    
end






