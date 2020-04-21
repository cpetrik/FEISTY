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

fCTsF = NaN*ones(nparam,95*12);
fCTsP = fCTsF;
fCTsD = fCTsF;
fCTmF = fCTsF;
fCTmP = fCTsF;
fCTmD = fCTsF;
fCTlP = fCTsF;
fCTlD = fCTsF;

fCSsF = NaN*ones(48111,nparam);
fCSsP = fCSsF;
fCSsD = fCSsF;
fCSmF = fCSsF;
fCSmP = fCSsF;
fCSmD = fCSsF;
fCSlP = fCSsF;
fCSlD = fCSsF;

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
    netcdf_read_fore_fished_cmax_ens(fname,simname);

    load([fname '_Means_cmax_' simname '.mat'],'yr50',...
        'sf_tcmax','sp_tcmax','sd_tcmax',...
        'mf_tcmax','mp_tcmax','md_tcmax',...
        'lp_tcmax','ld_tcmax',...
        'sf_cmax50','sp_cmax50','sd_cmax50',...
        'mf_cmax50','mp_cmax50','md_cmax50',...
        'lp_cmax50','ld_cmax50');

    %% Time series
    fCTsF(j,:) = sf_tcmax;
    fCTsP(j,:) = sp_tcmax;
    fCTsD(j,:) = sd_tcmax;
    fCTmF(j,:) = mf_tcmax;
    fCTmP(j,:) = mp_tcmax;
    fCTmD(j,:) = md_tcmax;
    fCTlP(j,:) = lp_tcmax;
    fCTlD(j,:) = ld_tcmax;

    %% Last 50 years (2051-2100)
    fCSsF(:,j) = sf_cmax50;
    fCSsP(:,j) = sp_cmax50;
    fCSsD(:,j) = sd_cmax50;
    fCSmF(:,j) = mf_cmax50;
    fCSmP(:,j) = mp_cmax50;
    fCSmD(:,j) = md_cmax50;
    fCSlP(:,j) = lp_cmax50;
    fCSlD(:,j) = ld_cmax50;

    %% Maps
%     map_fore_cmax_ensem(simname,sf_cmax50,sp_cmax50,sd_cmax50,...
%         mf_cmax50,mp_cmax50,md_cmax50,lp_cmax50,ld_cmax50,pp);

end
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Forecast_All_fish03_ensem6_pdiffMaxMin1951_cmaxs.mat'],...
    'fCTsF','fCTsP','fCTsD','fCTmF','fCTmP','fCTmD','fCTlP','fCTlD',...
    'fCSsF','fCSsP','fCSsD','fCSmF','fCSmP','fCSmD','fCSlP','fCSlD',...
    'fnms','snms')
