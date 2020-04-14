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

%% grid
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

load([cpath 'lme_mask_esm2m.mat']);
tlme = lme_mask_esm2m';
lme_grid = tlme(ID);

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
params = red_params;

nparam = length(params);

%time series of catch summed each year
Hlme_tsc_mf = NaN*ones(nparam+1,145,66);
Hlme_tsc_mp = Hlme_tsc_mf;
Hlme_tsc_md = Hlme_tsc_mf;
Hlme_tsc_lp = Hlme_tsc_mf;
Hlme_tsc_ld = Hlme_tsc_mf;

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
    load([fname '_Means_' simname '.mat'],...
        'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc');
    
    %% Time series of catch (g per km2 per year) in each LME
    for L=1:66
        lid = find(lme_grid==L);
        Hlme_tsc_mf(j,:,L) = nansum(mf_tyc(lid,:)); 
        Hlme_tsc_mp(j,:,L) = nansum(mp_tyc(lid,:)); 
        Hlme_tsc_md(j,:,L) = nansum(md_tyc(lid,:)); 
        Hlme_tsc_lp(j,:,L) = nansum(lp_tyc(lid,:)); 
        Hlme_tsc_ld(j,:,L) = nansum(ld_tyc(lid,:)); 
    end
    
end

%% Add baseline param set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
load([fpath 'Means_Historic_' harv '_' cfile '.mat'],...
    'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc');

for L=1:66
        lid = find(lme_grid==L);
        Hlme_tsc_mf(nparam+1,:,L) = nansum(mf_tyc(lid,:)); 
        Hlme_tsc_mp(nparam+1,:,L) = nansum(mp_tyc(lid,:)); 
        Hlme_tsc_md(nparam+1,:,L) = nansum(md_tyc(lid,:)); 
        Hlme_tsc_lp(nparam+1,:,L) = nansum(lp_tyc(lid,:)); 
        Hlme_tsc_ld(nparam+1,:,L) = nansum(ld_tyc(lid,:)); 
end

%%
units_catch = 'g_km2_yr';

epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'LME_ts_yield_Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'Hlme_tsc_mf','Hlme_tsc_mp','Hlme_tsc_md','Hlme_tsc_lp','Hlme_tsc_ld',...
    'units_catch')

