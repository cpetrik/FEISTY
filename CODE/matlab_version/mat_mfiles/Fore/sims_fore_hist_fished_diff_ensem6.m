% Which sims have biggest changes
% 6 most sens parameters

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Prod_diff_50yr_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat'])
load([epath 'simnames_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat']);

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt2_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
params = red_params;


%% which sims have biggest changes
big = NaN*ones(7,2);
for i=1:7
    id = find(fbar(:,i)==min(fbar(:,i)));
    big(i,1) = id;
    big(i,2) = min(fbar(:,i));
end

%% tables
Stab = array2table(big,'VariableNames',{'sim','decr'},...
    'RowNames',names);
writetable(Stab,[epath 'Hist_Fore_All_fish03_ensem_max_pdiff_sims.csv'],...
    'Delimiter',',','WriteRowNames',true)

big2 = big;
big2(:,3:8) = params(big(:,1),:);
tab = array2table(big2,'VariableNames',{'sim','decr','assim','bM','bE','aM','aE','kt'},...
    'RowNames',names);
writetable(tab,[epath 'Hist_Fore_All_fish03_ensem6_mid_kt2_max_pdiff_sims_params.csv'],...
    'Delimiter',',','WriteRowNames',true)







