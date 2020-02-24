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
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

area = AREA_OCN(ID);
area_mat = repmat(area,1,length(y));

pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Fore_param6_mid_best/'];
if (~isfolder(pp))
    mkdir(pp)
end

%% LTL info for TEeffs
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([gpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzloss_mean_fore',...
    'lzloss_mean_fore','det_mean_fore','npp_mean_fore','mzprod_mean_fore',...
    'lzprod_mean_fore'); 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mmz_loss = mzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_loss = lzloss_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mmz_prod = mzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_prod = lzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mdet = det_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mnpp = npp_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([nfile 'LHS_param6_mid6_kt2_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
%     'red_params');
% params = red_params;
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
params = red_params;
nparam = length(params);

fTsF = NaN*ones(nparam,95*12);
fTsP = fTsF;
fTsD = fTsF;
fTmF = fTsF;
fTmP = fTsF;
fTmD = fTsF;
fTlP = fTsF;
fTlD = fTsF;
fTB  = fTsF;

fSsF = NaN*ones(48111,nparam);
fSsP = fSsF;
fSsD = fSsF;
fSmF = fSsF;
fSmP = fSsF;
fSmD = fSsF;
fSlP = fSsF;
fSlD = fSsF;
fSB  = fSsF;

Flme_mcatch = NaN*ones(66,5,nparam);
Flme_mbio   = NaN*ones(66,9,nparam);

fTEeffM = NaN*ones(48111,nparam);
fTEeff_ATL = fTEeffM;
fTEeff_LTL = fTEeffM;
fTEeff_HTL = fTEeffM;

ftTEeff_LTL = NaN*ones(1,95/5);
ftTEeffM = NaN*ones(nparam,95/5);
ftTEeff_ATL = ftTEeffM;
ftTEeff_HTL = ftTEeffM;

ftTE_ATL = ftTEeffM;
ftTE_HTL = ftTEeffM;

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
    if (j>=27)
        netcdf_ts_areaint_fore_fished_bio_prod_ens(fname,simname,area_mat);
    end
    
    load([fname '_Means_' simname '.mat'],'yr50',...
        'sf_tamean','sp_tamean','sd_tamean',...
        'mf_tamean','mp_tamean','md_tamean',...
        'lp_tamean','ld_tamean','b_tamean');
    
    %% Time series
    fTsF(j,:) = sf_tamean;
    fTsP(j,:) = sp_tamean;
    fTsD(j,:) = sd_tamean;
    fTmF(j,:) = mf_tamean;
    fTmP(j,:) = mp_tamean;
    fTmD(j,:) = md_tamean;
    fTlP(j,:) = lp_tamean;
    fTlD(j,:) = ld_tamean;
    fTB(j,:)  = b_tamean;
    
    
end
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'fTsF','-append')
