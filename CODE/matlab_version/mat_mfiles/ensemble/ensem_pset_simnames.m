% Make list of ensemble names

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

%%
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param5_mid5_bestAIC_params_multFup_neg_multPneg.mat'],'params');
nparam = length(params);

pnames = cell(length(params),1);
snames = cell(length(params),1);

%%
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params5(pset)
    
    %! Make core parameters/constants (global)
    const_params5()
    
    %! Create a directory for output
    simname = sub_fname_ens(frate);
    
    cfile = ['/Volumes/GFDL/NC/Matlab_new_size/' simname];
    
    snames{j}=simname;
    pnames{j}=cfile;
    
end
%%
epath = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
save([epath 'simnames_ensem5_mid5_bestAIC_multFup_multPneg.mat'],...
    'pnames','snames')
