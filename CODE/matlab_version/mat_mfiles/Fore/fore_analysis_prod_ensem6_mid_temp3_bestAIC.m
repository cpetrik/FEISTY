% FEISTY Forecast runs of best parameter sets
% varying 6 most sensitive params

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
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Fore_param6_mid_best/'];
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
params = red_params;
nparam = length(params);

fSPsF = NaN*ones(48111,nparam);
fSPsP = fSPsF;
fSPsD = fSPsF;
fSPmF = fSPsF;
fSPmP = fSPsF;
fSPmD = fSPsF;
fSPlP = fSPsF;
fSPlD = fSPsF;
fSPB  = fSPsF;

%%
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params6(pset)
    
    %! Make core parameters/constants (global)
    const_params6()
    
    %! Create a directory for output
    [fname,simname] = sub_fname_fore_ens(frate);
    
    cfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' simname];

    %% Prod results for TEs
    load([fname '_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50');
    load([fname '_Means_' simname '.mat'],'yr50',...
        'b_mean50');
    
    % Last 50 years (2051-2100)
    fSPsF(:,j) = sf_prod50;
    fSPsP(:,j) = sp_prod50;
    fSPsD(:,j) = sd_prod50;
    fSPmF(:,j) = mf_prod50;
    fSPmP(:,j) = mp_prod50;
    fSPmD(:,j) = md_prod50;
    fSPlP(:,j) = lp_prod50;
    fSPlD(:,j) = ld_prod50;
    fSPB(:,j)  = b_mean50;
    
end
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
lpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'fSPsF','fSPsP','fSPsD','fSPmF','fSPmP','fSPmD','fSPlP','fSPlD','fSPB',...
    '-append')
save([lpath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'fSPsF','fSPsP','fSPsD','fSPmF','fSPmP','fSPmD','fSPlP','fSPlD','fSPB',...
    '-append')
