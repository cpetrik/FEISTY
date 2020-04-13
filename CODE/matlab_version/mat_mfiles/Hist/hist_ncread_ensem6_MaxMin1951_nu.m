% FEISTY Historic runs of parameter sets
% that produced the max and min changes in production

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

gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
load([nfile 'Hist_Fore_All_fish03_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],'pstats');

params = red_params;
psets = unique(pstats(:,1));
nparam = length(psets);

hNTsF = NaN*ones(nparam,145*12);
hNTsP = hNTsF;
hNTsD = hNTsF;
hNTmF = hNTsF;
hNTmP = hNTsF;
hNTmD = hNTsF;
hNTlP = hNTsF;
hNTlD = hNTsF;

hNSsF = NaN*ones(48111,nparam);
hNSsP = hNSsF;
hNSsD = hNSsF;
hNSmF = hNSsF;
hNSmP = hNSsF;
hNSmD = hNSsF;
hNSlP = hNSsF;
hNSlD = hNSsF;

%%
for j = 1:nparam
    %! Change individual parameters
    k = psets(j);
    
    %! Change individual parameters
    pset = params(k,:);
    set_params6(pset)
    
    %! Make core parameters/constants (global)
    const_params6()
    
    %! Create a directory for output
    [fname,simname] = sub_fname_hist_ens(frate);
    
    cfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' simname];
    
    %% Last 50 year means
%     netcdf_read_hist_fished_nu_ens(fname,simname)
    
    load([fname '_Means_nu_' simname '.mat'],'y','yr50',...
        'sf_tnu','sp_tnu','sd_tnu',...
        'mf_tnu','mp_tnu','md_tnu',...
        'lp_tnu','ld_tnu',...
        'sf_nu50','sp_nu50','sd_nu50',...
        'mf_nu50','mp_nu50','md_nu50',...
        'lp_nu50','ld_nu50');
    
    
    %% Time series
    hNTsF(j,:) = sf_tnu;
    hNTsP(j,:) = sp_tnu;
    hNTsD(j,:) = sd_tnu;
    hNTmF(j,:) = mf_tnu;
    hNTmP(j,:) = mp_tnu;
    hNTmD(j,:) = md_tnu;
    hNTlP(j,:) = lp_tnu;
    hNTlD(j,:) = ld_tnu;
    
    %% Last 50 years (1951-2000)
    hNSsF(:,j) = sf_nu50;
    hNSsP(:,j) = sp_nu50;
    hNSsD(:,j) = sd_nu50;
    hNSmF(:,j) = mf_nu50;
    hNSmP(:,j) = mp_nu50;
    hNSmD(:,j) = md_nu50;
    hNSlP(:,j) = lp_nu50;
    hNSlD(:,j) = ld_nu50;
    
    %% Maps
%     map_hist_nu_ensem(simname,sf_nu50,sp_nu50,sd_nu50,...
%         mf_nu50,mp_nu50,md_nu50,lp_nu50,ld_nu50,pp);
    
end
%%
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Historic_All_fish03_ensem6_MaxMin1951_nus.mat'],...
    'hNTsF','hNTsP','hNTsD','hNTmF','hNTmP','hNTmD','hNTlP','hNTlD',...
    'hNSsF','hNSsP','hNSsD','hNSmF','hNSmP','hNSmD','hNSlP','hNSlD')

