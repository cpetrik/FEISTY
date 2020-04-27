% FEISTY Forecast runs of best parameter sets
% varying 5 most sensitive params, same k's
% calc ts of area-int fractions

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
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat','AREA_OCN');
grid = csvread([cpath 'grid_csv.csv']);
area = AREA_OCN(grid(:,1));
area_mat = repmat(area,1,95*12);

pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Fore_param6_mid_best/'];

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_samek_bestAIC_params_Fupneg_mult8_Pneg2_mult3.mat'],...
    'params');
nparam = length(params);

faTsF = NaN*ones(nparam,95*12);
faTsP = faTsF;
faTsD = faTsF;
faTmF = faTsF;
faTmP = faTsF;
faTmD = faTsF;
faTlP = faTsF;
faTlD = faTsF;
faTB  = faTsF;

fatPD   = NaN*ones(nparam,95*12);
fatFD   = fatPD;
fatPelD = fatPD;
fatDPel = fatPD;
fatDP   = fatPD;
fatDF   = fatPD;
fatFP   = fatPD;
fatPF   = fatPD;
fatF    = fatPD;
fatP    = fatPD;
fatPel  = fatPD;

%%
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params6_samek(pset)
    
    %! Make core parameters/constants (global)
    const_params6_samek()
    
    %! Create a directory for output
    [fname,simname] = sub_fname_fore_ens_samek(frate);
    
    cfile = ['/Volumes/FEISTY/NC/Matlab_new_size/' simname];
    
    %% timeseries with area
    netcdf_ts_areaint_fore_fished_bio_prod_ens(fname,simname,area_mat);
    
    load([fname '_Means_' simname '.mat'],...
        'sf_tamean','sp_tamean','sd_tamean',...
        'mf_tamean','mp_tamean','md_tamean',...
        'lp_tamean','ld_tamean','b_tamean',...
    'tPD','tFD','tPelD','tDPel','tDP','tDF','tFP','tPF','tF','tP','tPel');
    
    %% Time series
    faTsF(j,:) = sf_tamean;
    faTsP(j,:) = sp_tamean;
    faTsD(j,:) = sd_tamean;
    faTmF(j,:) = mf_tamean;
    faTmP(j,:) = mp_tamean;
    faTmD(j,:) = md_tamean;
    faTlP(j,:) = lp_tamean;
    faTlD(j,:) = ld_tamean;
    faTB(j,:)  = b_tamean;
    
    fatPD(j,:) = tPD;
    fatFD(j,:) = tFD;
    fatPelD(j,:) = tPelD;
    fatDPel(j,:) = tDPel;
    fatDP(j,:) = tDP;
    fatDF(j,:) = tDF;
    fatFP(j,:) = tFP;
    fatPF(j,:) = tPF;
    fatF(j,:)  = tF;
    fatP(j,:)  = tP;
    fatPel(j,:)= tPel;
    
    
end
lpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

save([nfile 'Forecast_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'faTsF','faTsP','faTsD','faTmF','faTmP','faTmD','faTB','faTlP','faTlD',...
    'fatPD','fatFD','fatPelD','fatDPel','fatDP','fatDF','fatFP','fatPF',...
    'fatF','fatP','fatPel','-append')
save([lpath 'Forecast_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'faTsF','faTsP','faTsD','faTmF','faTmP','faTmD','faTB','faTlP','faTlD',...
    'fatPD','fatFD','fatPelD','fatDPel','fatDP','fatDF','fatFP','fatPF',...
    'fatF','fatP','fatPel','-append')
