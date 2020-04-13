% Table of pure change and percent change of
% production and fisheries yield

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
efile = 'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',cfile,'/full_runs/'];

%% Ensemble parameter sets
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
epath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];

%% take mean and error bars
load([epath 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])

ppmean = mean(fbar*100);
ppstd = std(fbar*100);

dpmean = mean(dbar);
dpstd = std(dbar);

prod_nm = names;

clear fbar dbar names

%%
load([epath 'Hist_Fore_All_fish03_ensem6_mid_kt3_yield_diff_pdiff.mat'],...
    'fbar','dbar','names')
% load([epath2 'Hist_Fore_All_fish03_ensem6_mid_kt3_yield_diff_pdiff.mat'],...
%     'fbar','dbar','names')

pymean = mean(fbar*100);
pystd = std(fbar*100);

dymean = mean(dbar*1e-6);
dystd = std(dbar*1e-6);

yield_nm = names;

clear fbar dbar names

%% tables
pstats(:,1) = dpmean';
pstats(:,2) = dpstd';
pstats(:,3) = ppmean';
pstats(:,4) = ppstd';

ystats(:,1) = dymean';
ystats(:,2) = dystd';
ystats(:,3) = pymean';
ystats(:,4) = pystd';

Ptab = array2table(pstats,'VariableNames',{'mean prod diff','std prod diff',...
    'mean prod pdiff','std prod pdiff'},'RowNames',prod_nm);
Ytab = array2table(ystats,'VariableNames',{'mean yield diff','std yield diff',...
    'mean yield pdiff','std yield pdiff'},'RowNames',yield_nm);

%%
writetable(Ptab,[epath 'Hist_Fore_All_fish03_ensem6_mid_kt3_prod_diff_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ptab,[epath2 'Hist_Fore_All_fish03_ensem6_mid_kt3_prod_diff_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)

writetable(Ytab,[epath 'Hist_Fore_All_fish03_ensem6_mid_kt3_yield_diff_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ytab,[epath2 'Hist_Fore_All_fish03_ensem6_mid_kt3_yield_diff_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)

