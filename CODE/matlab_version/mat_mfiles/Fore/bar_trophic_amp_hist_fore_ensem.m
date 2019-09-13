% Bar plot of trophic amplification

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

epath = '/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Prod_diff_50yr_ensem5_mid5_bestAIC_multFup_multPneg.mat'])

%% take mean and error bars
fmean = mean(fbar);
fstd = std(fbar);

%% bar graphs
figure(1)
bar(100*pbar)
set(gca,'XTickLabel',pnames)

%%
figure(2)
bar(100*fmean); hold on;
er = errorbar(1:7,100*fmean,100*fstd);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca,'XTickLabel',names)
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_all_types_ensem.png'])

%% size
sbar = pbar(1:2);
sbar(3:4) = fmean(1:2);
sstd = fstd(1:2);

figure(3)
bar(sbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),100*sstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','Mesozoo','M Fishes','L Fishes'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_ensem.png'])


%% type
tbar = fmean(3:5);
tstd = fstd(3:5);

figure(4)
bar(tbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(1:3,100*tbar,100*tstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_type_ensem.png'])

%% tables
tstats(:,1) = fmean';
tstats(:,2) = fstd';
Stab = array2table(tstats,'VariableNames',{'mean','std'},...
    'RowNames',names);
writetable(Stab,[epath 'Hist_Fore_All_fish03_ensem_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)







