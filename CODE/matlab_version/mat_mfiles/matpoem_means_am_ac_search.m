% POEM output at all locations
% Only last 12 months of 150 years saved 

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
%cfile = 'Dc_enc50-b210_m4-b175-k060_c50-b250_D075_J075_A050_Sm025_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';

%! Make core parameters/constants (global)
% make_parameters_mid()
% A = 0.5;
% bpow = 0.175;
% bcmx = 0.25;

aep = 10:10:100;
acp = 10:10:100;

SF = NaN*ones(length(aep),length(acp),16);
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
fSF = SF;
fSP = SF;
fSD = SF;
fMF = SF;
fMP = SF;
fMD = SF;
fLP = SF;
fLD = SF;
gSF = SF;
gSP = SF;
gSD = SF;
gMF = SF;
gMP = SF;
gMD = SF;
gLP = SF;
gLD = SF;

%%
for k=1:length(acp)
    for j=1:length(aep)   
        h = acp(k);         % coeff on Cmax
        gam = aep(j);       % coeff on search area

    
%     sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
%         '/param_sens/Climatol_All_fish03_enc',num2str(gam),'_c',num2str(h),'_locs.mat'];
    sfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,...
        '/param_sens/Climatol_All_fish03_enc-a',num2str(gam),'_c-a',num2str(h),'_locs.mat'];
    %Climatol_All_fish03_enc-a100_c-a100_locs
    load(sfile);
    
    % Last year
    [nt,nv,id] = size(S_Cobalt);
    time=1:nt;
    lyr=time((end-12+1):end);
    %biomass
    sp_mean=squeeze(mean(S_Sml_p(lyr,1,:),1));
    sf_mean=squeeze(mean(S_Sml_f(lyr,1,:),1));
    sd_mean=squeeze(mean(S_Sml_d(lyr,1,:),1));
    mp_mean=squeeze(mean(S_Med_p(lyr,1,:),1));
    mf_mean=squeeze(mean(S_Med_f(lyr,1,:),1));
    md_mean=squeeze(mean(S_Med_d(lyr,1,:),1));
    lp_mean=squeeze(mean(S_Lrg_p(lyr,1,:),1));
    ld_mean=squeeze(mean(S_Lrg_d(lyr,1,:),1));
    b_mean=squeeze(mean(S_Cobalt(lyr,1,:),1));
    
    %catch
    mf_my=squeeze(mean(S_Med_f(lyr,25,:),1));
    mp_my=squeeze(mean(S_Med_p(lyr,25,:),1));
    md_my=squeeze(mean(S_Med_d(lyr,25,:),1));
    lp_my=squeeze(mean(S_Lrg_p(lyr,25,:),1));
    ld_my=squeeze(mean(S_Lrg_d(lyr,25,:),1));
    
    %gut fullness (= con/cmax)
    sp_gut=squeeze(mean(S_Sml_p(lyr,20,:),1));
    sf_gut=squeeze(mean(S_Sml_f(lyr,20,:),1));
    sd_gut=squeeze(mean(S_Sml_d(lyr,20,:),1));
    mp_gut=squeeze(mean(S_Med_p(lyr,20,:),1));
    mf_gut=squeeze(mean(S_Med_f(lyr,20,:),1));
    md_gut=squeeze(mean(S_Med_d(lyr,20,:),1));
    lp_gut=squeeze(mean(S_Lrg_p(lyr,20,:),1));
    ld_gut=squeeze(mean(S_Lrg_d(lyr,20,:),1));
    
    %GGE
    sp_gge=squeeze(mean(S_Sml_p(lyr,15,:)./S_Sml_p(lyr,14,:),1));
    sf_gge=squeeze(mean(S_Sml_f(lyr,15,:)./S_Sml_f(lyr,14,:),1));
    sd_gge=squeeze(mean(S_Sml_d(lyr,15,:)./S_Sml_d(lyr,14,:),1));
    mp_gge=squeeze(mean(S_Med_p(lyr,15,:)./S_Med_p(lyr,14,:),1));
    mf_gge=squeeze(mean(S_Med_f(lyr,15,:)./S_Med_f(lyr,14,:),1));
    md_gge=squeeze(mean(S_Med_d(lyr,15,:)./S_Med_d(lyr,14,:),1));
    lp_gge=squeeze(mean(S_Lrg_p(lyr,15,:)./S_Lrg_p(lyr,14,:),1));
    ld_gge=squeeze(mean(S_Lrg_d(lyr,15,:)./S_Lrg_d(lyr,14,:),1));
    
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
    
    fSF(j,k,:) = sf_gut;
    fSP(j,k,:) = sp_gut;
    fSD(j,k,:) = sd_gut;
    fMF(j,k,:) = mf_gut;
    fMP(j,k,:) = mp_gut;
    fMD(j,k,:) = md_gut;
    fLP(j,k,:) = lp_gut;
    fLD(j,k,:) = ld_gut;
    
    gSF(j,k,:) = sf_gge;
    gSP(j,k,:) = sp_gge;
    gSD(j,k,:) = sd_gge;
    gMF(j,k,:) = mf_gge;
    gMP(j,k,:) = mp_gge;
    gMD(j,k,:) = md_gge;
    gLP(j,k,:) = lp_gge;
    gLD(j,k,:) = ld_gge;
    
    end
end

%%
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
save([nfile 'Locs_Climatol_All_fish03_means_aenc_acmax_search.mat'])






