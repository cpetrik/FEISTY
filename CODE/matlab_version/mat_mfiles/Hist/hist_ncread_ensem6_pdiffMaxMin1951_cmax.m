% FEISTY Historic runs of parameter sets
% that produced the max and min percent changes in production

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
load([nfile 'Hist_Fore_All_fish03_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],'pdstat');

params = red_params;
psets = unique(pdstat(:,1));
nparam = length(psets);

hCTsF = NaN*ones(nparam,145*12);
hCTsP = hCTsF;
hCTsD = hCTsF;
hCTmF = hCTsF;
hCTmP = hCTsF;
hCTmD = hCTsF;
hCTlP = hCTsF;
hCTlD = hCTsF;

hCSsF = NaN*ones(48111,nparam);
hCSsP = hCSsF;
hCSsD = hCSsF;
hCSmF = hCSsF;
hCSmP = hCSsF;
hCSmD = hCSsF;
hCSlP = hCSsF;
hCSlD = hCSsF;

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
    netcdf_read_hist_fished_cmax_ens(fname,simname)
    
    load([fname '_Means_cmax_' simname '.mat'],'y','yr50',...
        'sf_tcmax','sp_tcmax','sd_tcmax',...
        'mf_tcmax','mp_tcmax','md_tcmax',...
        'lp_tcmax','ld_tcmax',...
        'sf_cmax50','sp_cmax50','sd_cmax50',...
        'mf_cmax50','mp_cmax50','md_cmax50',...
        'lp_cmax50','ld_cmax50');
    
    
    %% Time series
    hCTsF(j,:) = sf_tcmax;
    hCTsP(j,:) = sp_tcmax;
    hCTsD(j,:) = sd_tcmax;
    hCTmF(j,:) = mf_tcmax;
    hCTmP(j,:) = mp_tcmax;
    hCTmD(j,:) = md_tcmax;
    hCTlP(j,:) = lp_tcmax;
    hCTlD(j,:) = ld_tcmax;
    
    %% Last 50 years (1951-2000)
    hCSsF(:,j) = sf_cmax50;
    hCSsP(:,j) = sp_cmax50;
    hCSsD(:,j) = sd_cmax50;
    hCSmF(:,j) = mf_cmax50;
    hCSmP(:,j) = mp_cmax50;
    hCSmD(:,j) = md_cmax50;
    hCSlP(:,j) = lp_cmax50;
    hCSlD(:,j) = ld_cmax50;
    
    %% Maps
%     map_hist_cmax_ensem(simname,sf_cmax50,sp_cmax50,sd_cmax50,...
%         mf_cmax50,mp_cmax50,md_cmax50,lp_cmax50,ld_cmax50,pp);
    
end
%%
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Historic_All_fish03_ensem6_pdiffMaxMin1951_cmaxs.mat'],...
    'hCTsF','hCTsP','hCTsD','hCTmF','hCTmP','hCTmD','hCTlP','hCTlD',...
    'hCSsF','hCSsP','hCSsD','hCSmF','hCSmP','hCSmD','hCSlP','hCSlD')

