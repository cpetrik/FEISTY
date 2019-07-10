% FEISTY Historic runs of best parameter sets 
% varying 5 most sensitive params

clear all
close all

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
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
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' ...
    'param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    'Hist_param5_mid5_best/'];

%%
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param5_mid5_bestAIC_params_multFup_neg_multPneg.mat'],'params');
nparam = length(params);

hTsF = NaN*ones(nparam,95);
hTsP = fTsF;
hTsD = fTsF;
hTmF = fTsF;
hTmP = fTsF;
hTmD = fTsF;
hTlP = fTsF;
hTlD = fTsF;
hTB  = fTsF;

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

for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params5(pset)
    
    %! Make core parameters/constants (global)
    const_params5()
    
    %! Create a directory for output
    [fname,simname] = sub_fname_hist_ens(frate);
    
    cfile = ['/Volumes/GFDL/NC/Matlab_new_size/' simname];
    
    %% Last year means
    netcdf_read_hist_fished_bio_ens(fname,simname)
    
    %% Time series
    hTsF(j,:) = sf_tmean;
    hTsP(j,:) = sp_tmean;
    hTsD(j,:) = sd_tmean;
    hTmF(j,:) = mf_tmean;
    hTmP(j,:) = mp_tmean;
    hTmD(j,:) = md_tmean;
    hTlP(j,:) = lp_tmean;
    hTlD(j,:) = ld_tmean;
    hTB(j,:)  = b_tmean;
    
    %% Last 50 years (2051-2100)
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
    %     vis_hist_ensem(sname,sf_mean50,sp_mean50,sd_mean50,...
    %     mf_mean50,mp_mean50,md_mean50,b_mean50,lp_mean50,ld_mean50,pp);
    
    %% LME
%     [lme_mcatch,lme_mbio,lme_area] = lme_hist_ensem(sf_mean,sp_mean,sd_mean,...
%         mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,lp_my,ld_my);
%     Hlme_mcatch(:,:,j) = lme_mcatch;
%     Hlme_mbio(:,:,j)   = lme_mbio;

    %% netcdf read prod results for TEs
    
    %% TE Effs
%     [TEeffM,TEeff_ATL,TEeff_LTLd,TEeff_HTLd] = hist_fished_effTEs_useDet_ensem(BE,mnpp,mdet,mmz_loss,mlz_loss,...
%     mf_prod50,mp_prod50,md_prod50,lp_prod50,ld_prod50)
%     fTEeffM(:,j) = TEeffM;
%     fTEeff_ATL(:,j) = TEeff_ATL;
%     fTEeff_LTL(:,j) = TEeff_LTLd;
%     fTEeff_HTL(:,j) = TEeff_HTLd;
    
end
save('hTsF','hTsP','hTsD','hTmF','hTmP','hTmD','hTlP','hTlD','hTB',...
    'hSsF','hSsP','hSsD','hSmF','hSmP','hSmD','hSlP','hSlD','hSB',...
    'Hlme_mcatch','Hlme_mbio','lme_area',...
    'hTEeffM','hTEeff_ATL','hTEeff_LTL','hTEeff_HTL')

