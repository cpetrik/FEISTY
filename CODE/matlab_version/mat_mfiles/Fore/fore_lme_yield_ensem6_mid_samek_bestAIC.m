% FEISTY Forecast runs of best parameter sets
% varying 6 most sensitive params (added kt)
% total fisheries yield in each LMEs, time series

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
    'param_ensemble/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Hist_param6_mid_best/'];

%% grid
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

load([cpath 'lme_mask_esm2m.mat']);
tlme = lme_mask_esm2m';
lme_grid = tlme(ID);

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_samek_bestAIC_params_Fupneg_mult8_Pneg2_mult3.mat'],...
    'params');

nparam = length(params);

%time series of catch summed each year
Flme_tsc_mf = NaN*ones(nparam+1,95,66);
Flme_tsc_mp = Flme_tsc_mf;
Flme_tsc_md = Flme_tsc_mf;
Flme_tsc_lp = Flme_tsc_mf;
Flme_tsc_ld = Flme_tsc_mf;

%%
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params6_samek(pset)
    
    %! Make core parameters/constants (global)
    const_params6_samek()
    
    %! Create a directory for output
    [fname,simname] = sub_fname_fore_ens_samek(frate);
    
    %% Last 50 year means
    load([fname '_Means_' simname '.mat'],...
        'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc');
    
    %% Time series of catch (g per km2 per year) in each LME
    for L=1:66
        lid = find(lme_grid==L);
        Flme_tsc_mf(j,:,L) = nansum(mf_tyc(lid,:)); 
        Flme_tsc_mp(j,:,L) = nansum(mp_tyc(lid,:)); 
        Flme_tsc_md(j,:,L) = nansum(md_tyc(lid,:)); 
        Flme_tsc_lp(j,:,L) = nansum(lp_tyc(lid,:)); 
        Flme_tsc_ld(j,:,L) = nansum(ld_tyc(lid,:)); 
    end
    
end
%% Add baseline param set
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
load([fpath 'Means_fore_' harv '_' cfile '.mat'],...
    'mf_tyc','mp_tyc','md_tyc','lp_tyc','ld_tyc');

for L=1:66
        lid = find(lme_grid==L);
        Flme_tsc_mf(nparam+1,:,L) = nansum(mf_tyc(lid,:)); 
        Flme_tsc_mp(nparam+1,:,L) = nansum(mp_tyc(lid,:)); 
        Flme_tsc_md(nparam+1,:,L) = nansum(md_tyc(lid,:)); 
        Flme_tsc_lp(nparam+1,:,L) = nansum(lp_tyc(lid,:)); 
        Flme_tsc_ld(nparam+1,:,L) = nansum(ld_tyc(lid,:)); 
end
    
%%
units_catch = 'g_km2_yr';

lpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([lpath 'LME_ts_yield_Forecast_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'Flme_tsc_mf','Flme_tsc_mp','Flme_tsc_md','Flme_tsc_lp','Flme_tsc_ld',...
    'units_catch')
save([nfile 'LME_ts_yield_Forecast_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'Flme_tsc_mf','Flme_tsc_mp','Flme_tsc_md','Flme_tsc_lp','Flme_tsc_ld',...
    'units_catch')

