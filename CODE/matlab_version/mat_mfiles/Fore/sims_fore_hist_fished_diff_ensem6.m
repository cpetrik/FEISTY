% Which sims have biggest changes
% 6 most sens parameters

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([epath 'Prod_diff_50yr_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat'])
% load([epath 'simnames_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat']);

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])
load([nfile 'simnames_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat']);
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
params = red_params;


%% which sims have biggest changes
%fbar = % diffs
%col = {'M','L','F','P','D','Pel','All','B'};
big = NaN*ones(8,2);
for i=1:8
    id = find(abs(fbar(:,i))==max(abs(fbar(:,i))));
    big(i,1) = id;
    big(i,2) = fbar(id,i);
end

sml = NaN*ones(8,2);
for i=1:8
    id2 = find(abs(fbar(:,i))==min(abs(fbar(:,i))));
    sml(i,1) = id2;
    sml(i,2) = fbar(id2,i);
end

%% tables
Stab = array2table(big,'VariableNames',{'sim','pdiff'},...
    'RowNames',names);
writetable(Stab,[nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_max_pdiff_sims.csv'],...
    'Delimiter',',','WriteRowNames',true)

big2 = big;
big2(:,3:8) = params(big(:,1),:);
tab = array2table(big2,'VariableNames',{'sim','pdiff','assim','bM','bE','aM','aE','kt'},...
    'RowNames',names);
writetable(tab,[nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_max_pdiff_sims_params.csv'],...
    'Delimiter',',','WriteRowNames',true)


Mtab = array2table(sml,'VariableNames',{'sim','pdiff'},...
    'RowNames',names);
writetable(Mtab,[nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_min_pdiff_sims.csv'],...
    'Delimiter',',','WriteRowNames',true)

sml2 = sml;
sml2(:,3:8) = params(sml(:,1),:);
tsm = array2table(sml2,'VariableNames',{'sim','pdiff','assim','bM','bE','aM','aE','kt'},...
    'RowNames',names);
writetable(tsm,[nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_min_pdiff_sims_params.csv'],...
    'Delimiter',',','WriteRowNames',true)

save([nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_maxmin_pdiffs.mat'],...
    'big2','sml2','tab','tsm','names','pnames','snames')







