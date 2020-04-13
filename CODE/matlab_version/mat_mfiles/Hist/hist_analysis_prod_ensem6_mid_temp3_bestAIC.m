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

hSPsF = NaN*ones(48111,nparam);
hSPsP = hSPsF;
hSPsD = hSPsF;
hSPmF = hSPsF;
hSPmP = hSPsF;
hSPmD = hSPsF;
hSPlP = hSPsF;
hSPlD = hSPsF;
hSPB  = hSPsF;

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
    load([fname '_Means_' simname '.mat'],'y','yr50',...
        'b_mean50');
    
    load([fname '_Means_prod_' simname '.mat'],...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50');
    
    %% Last 50 years (1951-2000)
    hSPsF(:,j) = sf_prod50;
    hSPsP(:,j) = sp_prod50;
    hSPsD(:,j) = sd_prod50;
    hSPmF(:,j) = mf_prod50;
    hSPmP(:,j) = mp_prod50;
    hSPmD(:,j) = md_prod50;
    hSPlP(:,j) = lp_prod50;
    hSPlD(:,j) = ld_prod50;
    hSPB(:,j)  = b_mean50;
    
end
%%
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
lpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'hSPsF','hSPsP','hSPsD','hSPmF','hSPmP','hSPmD','hSPlP','hSPlD','hSPB',...
    '-append')
save([lpath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'hSPsF','hSPsP','hSPsD','hSPmF','hSPmP','hSPmD','hSPlP','hSPlD','hSPB',...
    '-append')

