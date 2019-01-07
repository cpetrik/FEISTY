% POEM output at all locations
% Only last 12 months of 150 years saved 

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
   
aparam = 0.5:0.1:1;
jparam = 0.5:0.1:1;

SF = NaN*ones(length(jparam),length(aparam),16);
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

for j=1:length(jparam)
    for k=1:length(aparam)
        J = jparam(j);
        A = aparam(k);
        tj = num2str(1000+int64(100*J));
        ta = num2str(1000+int64(100*A));

    
    sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/param_sens/Climatol_All_fish03_J',tj(2:end),'_A',ta(2:end),'_locs.mat'];
    load(sfile);
    
    % Last year
    [nt,nv,id] = size(S_Cobalt);
    time=1:nt;
    lyr=time((end-12+1):end);
    sp_mean=squeeze(mean(S_Sml_p(lyr,1,:),1));
    sf_mean=squeeze(mean(S_Sml_f(lyr,1,:),1));
    sd_mean=squeeze(mean(S_Sml_d(lyr,1,:),1));
    mp_mean=squeeze(mean(S_Med_p(lyr,1,:),1));
    mf_mean=squeeze(mean(S_Med_f(lyr,1,:),1));
    md_mean=squeeze(mean(S_Med_d(lyr,1,:),1));
    lp_mean=squeeze(mean(S_Lrg_p(lyr,1,:),1));
    ld_mean=squeeze(mean(S_Lrg_d(lyr,1,:),1));
    b_mean=squeeze(mean(S_Cobalt(lyr,1,:),1));
    
    mf_my=squeeze(mean(S_Med_f(lyr,25,:),1));
    mp_my=squeeze(mean(S_Med_p(lyr,25,:),1));
    md_my=squeeze(mean(S_Med_d(lyr,25,:),1));
    lp_my=squeeze(mean(S_Lrg_p(lyr,25,:),1));
    ld_my=squeeze(mean(S_Lrg_d(lyr,25,:),1));
    
    %
    save(sfile,...
        'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
        'lp_mean','ld_mean','time','lyr',...
        'mf_my','mp_my','md_my','lp_my','ld_my',...
        '-append');
    
    SF(j,k,:) = sf_mean;
    SP(j,k,:) = sp_mean;
    SD(j,k,:) = sd_mean;
    MF(j,k,:) = mf_mean;
    MP(j,k,:) = mp_mean;
    MD(j,k,:) = md_mean;
    LP(j,k,:) = lp_mean;
    LD(j,k,:) = ld_mean;
    BI(j,k,:) = b_mean;
    
    MFc(j,k,:) = mf_my;
    MPc(j,k,:) = mp_my;
    MDc(j,k,:) = md_my;
    LPc(j,k,:) = lp_my;
    LDc(j,k,:) = ld_my;
    
    end
end

%%
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
save([nfile 'Locs_Climatol_All_fish03_means_AJsearch.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc',...
    'aparam','jparam')






