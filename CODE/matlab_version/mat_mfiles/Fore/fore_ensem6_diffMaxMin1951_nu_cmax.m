% FEISTY Forecast runs of parameter sets
% that produced the max and min changes in production
% plot nu scaled by cmax, both from FEISTY output

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

grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

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

load([epath 'Historic_All_fish03_ensem6_diffMaxMin1951_cmaxs.mat'],...
    'hCSsF','hCSsP','hCSsD','hCSmF','hCSmP','hCSmD','hCSlP','hCSlD')
load([epath 'Forecast_All_fish03_ensem6_diffMaxMin1951_cmaxs.mat'],...
    'fCSsF','fCSsP','fCSsD','fCSmF','fCSmP','fCSmD','fCSlP','fCSlD');

%% Last 50 years (2051-2100) 
%scale nu by cmax
 
HsnuSF = hNSsF ./ hCSsF;
HsnuSP = hNSsP ./ hCSsP;
HsnuSD = hNSsD ./ hCSsD;
HsnuMF = hNSmF ./ hCSmF;
HsnuMP = hNSmP ./ hCSmP;
HsnuMD = hNSmD ./ hCSmD;
HsnuLP = hNSlP ./ hCSlP;
HsnuLD = hNSlD ./ hCSlD;

FsnuSF = fNSsF ./ fCSsF;
FsnuSP = fNSsP ./ fCSsP;
FsnuSD = fNSsD ./ fCSsD;
FsnuMF = fNSmF ./ fCSmF;
FsnuMP = fNSmP ./ fCSmP;
FsnuMD = fNSmD ./ fCSmD;
FsnuLP = fNSlP ./ fCSlP;
FsnuLD = fNSlD ./ fCSlD;

%% Maps
for j = 1:nparam
    sname = ['diff_' snms{j}];
    map_hist_scnu_ensem(sname,HsnuSF(:,j),HsnuSP(:,j),HsnuSD(:,j),...
        HsnuMF(:,j),HsnuMP(:,j),HsnuMD(:,j),HsnuLP(:,j),HsnuLD(:,j),hpp);
    
    map_fore_scnu_ensem(sname,FsnuSF(:,j),FsnuSP(:,j),FsnuSD(:,j),...
        FsnuMF(:,j),FsnuMP(:,j),FsnuMD(:,j),FsnuLP(:,j),FsnuLD(:,j),fpp);
end



