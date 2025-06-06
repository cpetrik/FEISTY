% FEISTY Historic runs of best parameter sets
% varying 6 most sensitive params (added kt)

clear all
close all

global GRD
global PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global dfrate frate

cfn=nan;
efn=nan;
mfn=nan;

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Hist_param6_mid_best/'];
if (~isfolder(pp))
    mkdir(pp)
end

%% load mnpp,mdet,mmz_loss,mlz_loss from Hist
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([gpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzloss_mean_hist',...
    'lzloss_mean_hist','det_mean_hist','npp_mean_hist','mzprod_mean_hist',...
    'lzprod_mean_hist'); 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mmz_loss = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_loss = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mmz_prod = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_prod = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mdet = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mnpp = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;


%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([nfile 'LHS_param6_mid6_kt2_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
%     'red_params');
% params = red_params;
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
% newp = find(red_params(:,6)==0.0955);
% params = red_params(newp,:);
params = red_params;

nparam = length(params);

hTsF = NaN*ones(nparam,145*12);
hTsP = hTsF;
hTsD = hTsF;
hTmF = hTsF;
hTmP = hTsF;
hTmD = hTsF;
hTlP = hTsF;
hTlD = hTsF;
hTB  = hTsF;

hSsF = NaN*ones(48111,nparam);
hSsP = hSsF;
hSsD = hSsF;
hSmF = hSsF;
hSmP = hSsF;
hSmD = hSsF;
hSlP = hSsF;
hSlD = hSsF;
hSB  = hSsF;

Hlme_mcatch = NaN*ones(66,5,nparam);
Hlme_mbio   = NaN*ones(66,9,nparam);

hTEeffM = NaN*ones(48111,nparam);
hTEeff_ATL = hTEeffM;
hTEeff_LTL = hTEeffM;
hTEeff_HTL = hTEeffM;

htTEeff_LTL = NaN*ones(1,145/5);
htTEeffM = NaN*ones(nparam,145/5);
htTEeff_ATL = htTEeffM;
htTEeff_HTL = htTEeffM;

htTE_ATL = htTEeffM;
htTE_HTL = htTEeffM;

%%
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params6(pset)
    
    %! Make core parameters/constants (global)
    const_params6()
    
    %! Create a directory for output
    [fname,simname] = sub_fname_hist_ens(frate);
    
    cfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' simname];
    
    %% Last 50 year means
%     if (j>=21)
%          netcdf_read_hist_fished_bio_ens(fname,simname)
%     end
    
    load([fname '_Means_' simname '.mat'],'y','yr50',...
        'sf_tmean','sp_tmean','sd_tmean',...
        'mf_tmean','mp_tmean','md_tmean',...
        'lp_tmean','ld_tmean','b_tmean',...
        'sf_mean50','sp_mean50','sd_mean50',...
        'mf_mean50','mp_mean50','md_mean50',...
        'lp_mean50','ld_mean50','b_mean50',...
        'mf_my50','mp_my50','md_my50','lp_my50','ld_my50');
    
    
    %% Time series
    hTsF(j,:) = sf_tmean; %j=28, y(1081), 1.9457e+35
    hTsP(j,:) = sp_tmean; %j=28, y(1081), 1.9457e+35
    hTsD(j,:) = sd_tmean; %j=28, y(1081), 1.9457e+35
    hTmF(j,:) = mf_tmean;
    hTmP(j,:) = mp_tmean;
    hTmD(j,:) = md_tmean;
    hTlP(j,:) = lp_tmean;
    hTlD(j,:) = ld_tmean;
    hTB(j,:)  = b_tmean;
    
    %% Last 50 years (1951-2000)
    hSsF(:,j) = sf_mean50;
    hSsP(:,j) = sp_mean50;
    hSsD(:,j) = sd_mean50;
    hSmF(:,j) = mf_mean50;
    hSmP(:,j) = mp_mean50;
    hSmD(:,j) = md_mean50;
    hSlP(:,j) = lp_mean50;
    hSlD(:,j) = ld_mean50;
    hSB(:,j)  = b_mean50;
    
    %% Maps
    vis_hist_ensem(simname,sf_mean50,sp_mean50,sd_mean50,...
        mf_mean50,mp_mean50,md_mean50,b_mean50,lp_mean50,ld_mean50,pp);
    
    %% LME
    [lme_mcatch,lme_mbio,lme_area] = lme_hist_ensem(sf_mean50,sp_mean50,...
        sd_mean50,mf_mean50,mp_mean50,md_mean50,b_mean50,lp_mean50,ld_mean50,...
        mf_my50,mp_my50,md_my50,lp_my50,ld_my50,fname,simname);
    Hlme_mcatch(:,:,j) = lme_mcatch;
    Hlme_mbio(:,:,j)   = lme_mbio;
    
    %% netcdf read prod results for TEs
%     if (j>=15)
%          netcdf_read_hist_fished_prod_ens(fname,simname)
%     end
    
    load([fname '_Means_prod_' simname '.mat'],...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50',...
    'mf_prod','mp_prod','md_prod',...
    'lp_prod','ld_prod');
    
    %% TE Effs
    % Mean of last 50 yrs, all locations
%     [TEeffM,TEeff_ATL,TEeff_LTLd,TEeff_HTLd] = ...
%         hist_fished_effTEs_useDet_ensem(bent_eff,mnpp,mdet,mmz_loss,mlz_loss,...
%         mf_prod50,mp_prod50,md_prod50,lp_prod50,ld_prod50,fname,simname);
    [TEeffM,TEeff_ATL,TEeff_LTLd,TEeff_HTLd] = ...
        hist_fished_effTEs_Det_Zprod_ensem(bent_eff,mnpp,mdet,mmz_prod,mlz_prod,...
        mf_prod50,mp_prod50,md_prod50,lp_prod50,ld_prod50,fname,simname);
    hTEeffM(:,j) = TEeffM;
    hTEeff_ATL(:,j) = TEeff_ATL;
    hTEeff_LTL(:,j) = TEeff_LTLd;
    hTEeff_HTL(:,j) = TEeff_HTLd;
    
    % 5 yr Means, all locations
%     [TEeffM,TEeff_A,TEeff_LTL,TEeff_H] = ...
%         hist_fished_effTEs_useDet_ensem(bent_eff,mnpp,mdet,mmz_loss,mlz_loss,...
%         mf_prod,mp_prod,md_prod,lp_prod,ld_prod,fname,simname);
    [TEeffM,TEeff_A,TEeff_LTL,TEeff_H] = ...
        hist_fished_effTEs_Det_Zprod_ensem(bent_eff,mnpp,mdet,mmz_prod,mlz_prod,...
        mf_prod,mp_prod,md_prod,lp_prod,ld_prod,fname,simname);
    htTEeffM(j,:)    = mean(TEeffM);
    htTEeff_LTL      = mean(TEeff_LTL);
    htTEeff_ATL(j,:) = mean(TEeff_A);
    htTEeff_HTL(j,:) = mean(TEeff_H);
    
    htTE_ATL(j,:) = mean(TEeff_A.^(1/4));
    htTE_HTL(j,:) = mean(TEeff_H.^(1/3));
    
end
%%
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'hTsF','hTsP','hTsD','hTmF','hTmP','hTmD','hTlP','hTlD','hTB',...
    'hSsF','hSsP','hSsD','hSmF','hSmP','hSmD','hSlP','hSlD','hSB',...
    'Hlme_mcatch','Hlme_mbio','lme_area',...
    'hTEeffM','hTEeff_ATL','hTEeff_LTL','hTEeff_HTL',...
    'htTEeffM','htTEeff_ATL','htTEeff_LTL','htTEeff_HTL',...
    'htTE_ATL','htTE_HTL')

