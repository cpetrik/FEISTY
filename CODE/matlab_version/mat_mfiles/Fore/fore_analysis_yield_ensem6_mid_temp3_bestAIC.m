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

%time series of catch summed each year
Flme_tsc_mf = NaN*ones(nparam,95);
Flme_tsc_mp = Flme_tsc_mf;
Flme_tsc_md = Flme_tsc_mf;
Flme_tsc_lp = Flme_tsc_mf;
Flme_tsc_ld = Flme_tsc_mf;

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
    
    %% Last 50 year means
    netcdf_read_fore_fished_yield_ens(fname,simname)
    
    load([fname '_Means_' simname '.mat'],...
        'mf_tsyc','mp_tsyc','md_tsyc','lp_tsyc','ld_tsyc');
    
    %% Time series of catch (g per km2 per year)
    Flme_tsc_mf(j,:) = mf_tsyc; 
    Flme_tsc_mp(j,:) = mp_tsyc;
    Flme_tsc_md(j,:) = md_tsyc;
    Flme_tsc_lp(j,:) = lp_tsyc;
    Flme_tsc_ld(j,:) = ld_tsyc;
    
end
%%
units_catch = 'g_km2_yr';

epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Flme_tsc_mf','Flme_tsc_mp','Flme_tsc_md','Flme_tsc_lp','Flme_tsc_ld',...
    'units_catch','-append')

