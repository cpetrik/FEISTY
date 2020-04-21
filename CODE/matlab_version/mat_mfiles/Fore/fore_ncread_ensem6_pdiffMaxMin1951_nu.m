% FEISTY Forecast runs of parameter sets
% that produced the max and min percent changes in production

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
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Fore_param6_mid_best/'];
if (~isfolder(pp))
    mkdir(pp)
end

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
load([nfile 'Hist_Fore_All_fish03_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],'pdstat');

params = red_params;
psets = unique(pdstat(:,1));
nparam = length(psets);

fNTsF = NaN*ones(nparam,95*12);
fNTsP = fNTsF;
fNTsD = fNTsF;
fNTmF = fNTsF;
fNTmP = fNTsF;
fNTmD = fNTsF;
fNTlP = fNTsF;
fNTlD = fNTsF;

fNSsF = NaN*ones(48111,nparam);
fNSsP = fNSsF;
fNSsD = fNSsF;
fNSmF = fNSsF;
fNSmP = fNSsF;
fNSmD = fNSsF;
fNSlP = fNSsF;
fNSlD = fNSsF;

fnms = cell(nparam,1);
snms = cell(nparam,1);

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
    [fname,simname] = sub_fname_fore_ens(frate);
    fnms{j} = fname;
    snms{j} = simname;

    cfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' simname];

    %% Last 50 year means
    netcdf_read_fore_fished_nu_ens(fname,simname);

    load([fname '_Means_nu_' simname '.mat'],'yr50',...
        'sf_tnu','sp_tnu','sd_tnu',...
        'mf_tnu','mp_tnu','md_tnu',...
        'lp_tnu','ld_tnu',...
        'sf_nu50','sp_nu50','sd_nu50',...
        'mf_nu50','mp_nu50','md_nu50',...
        'lp_nu50','ld_nu50');

    %% Time series
    fNTsF(j,:) = sf_tnu;
    fNTsP(j,:) = sp_tnu;
    fNTsD(j,:) = sd_tnu;
    fNTmF(j,:) = mf_tnu;
    fNTmP(j,:) = mp_tnu;
    fNTmD(j,:) = md_tnu;
    fNTlP(j,:) = lp_tnu;
    fNTlD(j,:) = ld_tnu;

    %% Last 50 years (2051-2100)
    fNSsF(:,j) = sf_nu50;
    fNSsP(:,j) = sp_nu50;
    fNSsD(:,j) = sd_nu50;
    fNSmF(:,j) = mf_nu50;
    fNSmP(:,j) = mp_nu50;
    fNSmD(:,j) = md_nu50;
    fNSlP(:,j) = lp_nu50;
    fNSlD(:,j) = ld_nu50;

    %% Maps
%     map_fore_nu_ensem(simname,sf_nu50,sp_nu50,sd_nu50,...
%         mf_nu50,mp_nu50,md_nu50,lp_nu50,ld_nu50,pp);

end
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Forecast_All_fish03_ensem6_pdiffMaxMin1951_nus.mat'],...
    'fNTsF','fNTsP','fNTsD','fNTmF','fNTmP','fNTmD','fNTlP','fNTlD',...
    'fNSsF','fNSsP','fNSsD','fNSmF','fNSmP','fNSmD','fNSlP','fNSlD',...
    'fnms','snms')
