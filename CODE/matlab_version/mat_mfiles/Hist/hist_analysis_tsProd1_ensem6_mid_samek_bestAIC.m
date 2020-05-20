% FEISTY Historic runs of best parameter sets
% varying 5 most sensitive params plus kt
% all k's equal
% compile all ensemble values of 1-yr prod means per grid cell
% then calc 1-yr mean TEeffs

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
    'param_ensemble/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Hist_param6_mid_best/'];
if (~isfolder(pp))
    mkdir(pp)
end

%% load mnpp,mdet,mmz_loss,mlz_loss from Hist
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([gpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzloss_mean_hist',...
    'lzloss_mean_hist','det_mean_hist','npp_mean_hist','mzprod_mean_hist',...
    'lzprod_mean_hist');

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mmz_loss = mzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_loss = lzloss_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mmz_prod = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mlz_prod = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mdet = det_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mnpp = npp_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;


%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_samek_bestAIC_params_Fupneg_mult8_Pneg2_mult3.mat'],...
    'params');

nparam = length(params);

ht1TEeff_LTL = NaN*ones(1,145);
ht1TEeff_M = NaN*ones(nparam,145);
ht1TEeff_ATL = ht1TEeff_M;
ht1TEeff_HTL = ht1TEeff_M;

ht1TE_ATL = ht1TEeff_M;
ht1TE_HTL = ht1TEeff_M;

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);
hLprod1 = NaN*ones(length(ID),145,nparam);

%%
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params6_samek(pset)

    %! Make core parameters/constants (global)
    const_params6_samek()

    %! Create a directory for output
    [fname,simname] = sub_fname_hist_ens_samek(frate);

    cfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' simname];

    %% Prod results for TEs
    %netcdf_read_hist_fished_ts_prod1_ens(fname,simname);

    load([fname '_Means_prod_' simname '.mat'],...
        'mf_prod1','mp_prod1','md_prod1',...
        'lp_prod1','ld_prod1');

    hLprod1(:,:,j) = lp_prod1 + ld_prod1;

    %% TE Effs
    % 1 yr Means, all locations
    [TEeffM,TEeff_A,TEeff_LTL,TEeff_H] = ...
        hist_fished_effTEs_Det_Zprod_ensem(bent_eff,mnpp,mdet,mmz_prod,...
        mlz_prod,mf_prod1,mp_prod1,md_prod1,lp_prod1,ld_prod1,fname,simname);
    ht1TEeff_M(j,:)   = mean(TEeffM);
    ht1TEeff_LTL      = mean(TEeff_LTL);
    ht1TEeff_ATL(j,:) = mean(TEeff_A);
    ht1TEeff_HTL(j,:) = mean(TEeff_H);

    ht1TE_ATL(j,:) = mean(TEeff_A.^(1/4));
    ht1TE_HTL(j,:) = mean(TEeff_H.^(1/3));

end
%%
lpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([nfile 'Historic_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'ht1TEeff_M','ht1TEeff_ATL','ht1TEeff_LTL','ht1TEeff_HTL',...
    'ht1TE_ATL','ht1TE_HTL','hLprod1','-append')

save([lpath 'Historic_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'ht1TEeff_M','ht1TEeff_ATL','ht1TEeff_LTL','ht1TEeff_HTL',...
    'ht1TE_ATL','ht1TE_HTL','hLprod1','-append')
