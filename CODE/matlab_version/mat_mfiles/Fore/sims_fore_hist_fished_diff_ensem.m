% Which sims have biggest changes

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

epath = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Prod_diff_50yr_ensem5_mid5_bestAIC_multFup_multPneg.mat'])

nfile = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'LHS_param5_mid5_bestAIC_params_multFup_neg_multPneg.mat'],'params');
load([nfile 'simnames_ensem5_mid5_bestAIC_multFup_multPneg.mat']);


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
big2(:,3:7) = params(big(:,1),:);
tab = array2table(big2,'VariableNames',{'sim','decr','assim','bM','bE','aM','aE'},...
    'RowNames',names);
writetable(tab,[epath 'Hist_Fore_All_fish03_ensem_max_pdiff_sims_params.csv'],...
    'Delimiter',',','WriteRowNames',true)







