% FEISTY Forecast runs of parameter sets
% that produced the max and min changes in production

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

hpp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Hist_param6_mid_best/'];
fpp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/'...
    'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'...
    '/full_runs/Fore_param6_mid_best/'];
if (~isfolder(hpp))
    mkdir(hpp)
end
if (~isfolder(fpp))
    mkdir(fpp)
end

%%
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'ptemp_mean_hist',...
    'ptemp_mean_fore','btemp_mean_hist','btemp_mean_fore');

grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

ptemp_hist = ptemp_mean_hist(ID);
ptemp_fore = ptemp_mean_fore(ID);
btemp_hist = btemp_mean_hist(ID);
btemp_fore = btemp_mean_fore(ID);

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
load([nfile 'Hist_Fore_All_fish03_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],'pstats');

params = red_params;
psets = unique(pstats(:,1));
nparam = length(psets);

%%
epath = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Historic_All_fish03_ensem6_MaxMin1951_nus.mat'],...
    'hNSsF','hNSsP','hNSsD','hNSmF','hNSmP','hNSmD','hNSlP','hNSlD')
load([epath 'Forecast_All_fish03_ensem6_MaxMin1951_nus.mat'],...
    'fNSsF','fNSsP','fNSsD','fNSmF','fNSmP','fNSmD','fNSlP','fNSlD',...
    'fnms','snms');

%% Last 50 years (2051-2100)
%scale nu by cmax
%calc bootleg cmax from temp and size
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

HScmax = (exp(0.063.*(ptemp_hist-10.0)) .* 20 .* M_s.^(-0.25)) ./365.0;
HMPcmax = (exp(0.063.*(ptemp_hist-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
HLPcmax = (exp(0.063.*(ptemp_hist-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;
HMBcmax = (exp(0.063.*(btemp_hist-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
HLBcmax = (exp(0.063.*(btemp_hist-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;

FScmax = (exp(0.063.*(ptemp_fore-10.0)) .* 20 .* M_s.^(-0.25)) ./365.0;
FMPcmax = (exp(0.063.*(ptemp_fore-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
FLPcmax = (exp(0.063.*(ptemp_fore-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;
FMBcmax = (exp(0.063.*(btemp_fore-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
FLBcmax = (exp(0.063.*(btemp_fore-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;

%% 
HsnuSF = hNSsF ./ repmat(HScmax,1,nparam);
HsnuSP = hNSsP ./ repmat(HScmax,1,nparam);
HsnuSD = hNSsD ./ repmat(HScmax,1,nparam);
HsnuMF = hNSmF ./ repmat(HMPcmax,1,nparam);
HsnuMP = hNSmP ./ repmat(HMPcmax,1,nparam);
HsnuMD = hNSmD ./ repmat(HMBcmax,1,nparam);
HsnuLP = hNSlP ./ repmat(HLPcmax,1,nparam);
HsnuLD = hNSlD ./ repmat(HLBcmax,1,nparam);

FsnuSF = fNSsF ./ repmat(FScmax,1,nparam);
FsnuSP = fNSsP ./ repmat(FScmax,1,nparam);
FsnuSD = fNSsD ./ repmat(FScmax,1,nparam);
FsnuMF = fNSmF ./ repmat(FMPcmax,1,nparam);
FsnuMP = fNSmP ./ repmat(FMPcmax,1,nparam);
FsnuMD = fNSmD ./ repmat(FMBcmax,1,nparam);
FsnuLP = fNSlP ./ repmat(FLPcmax,1,nparam);
FsnuLD = fNSlD ./ repmat(FLBcmax,1,nparam);

%% Maps
for j = 1:nparam
    map_hist_scnu_ensem(snms{j},HsnuSF(:,j),HsnuSP(:,j),HsnuSD(:,j),...
        HsnuMF(:,j),HsnuMP(:,j),HsnuMD(:,j),HsnuLP(:,j),HsnuLD(:,j),hpp);
    
    map_fore_scnu_ensem(snms{j},FsnuSF(:,j),FsnuSP(:,j),FsnuSD(:,j),...
        FsnuMF(:,j),FsnuMP(:,j),FsnuMD(:,j),FsnuLP(:,j),FsnuLD(:,j),fpp);
end



