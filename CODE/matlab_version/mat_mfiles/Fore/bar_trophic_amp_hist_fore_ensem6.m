% Bar plot of trophic amplification
% with std dev from ensemble sims
% benthic production = benthic biomass * detritus flux * benthic efficiency

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([epath 'Prod_diff_50yr_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat'])
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])

%% take mean and error bars
%means of global means, should they be means of all cells???
fmean = mean(fbar);
fstd = std(fbar);
fmax = max(fbar);
fmin = min(fbar);

%% bar graphs
figure(1)
bar(100*pbar,'k')
set(gca,'XTickLabel',pnames)

%% 1 st dev
figure(2)
b=bar(100*fmean); hold on;
b.FaceColor = [0 0.5 0.75];
er = errorbar(1:8,100*fmean,100*fstd);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca,'XTickLabel',names)
ylabel('Percent change')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_all_types_ensem_1std.png'])

%% 2 st dev
figure(3)
b=bar(100*fmean); hold on;
b.FaceColor = [0 0.5 0.75];
er = errorbar(1:8,100*fmean,200*fstd);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca,'XTickLabel',names)
ylabel('Percent change')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_all_types_ensem_2std.png'])

%% min max
figure(4)
b=bar(100*fmean); hold on;
b.FaceColor = [0 0.5 0.75];
er = errorbar(1:8,100*fmean,100*(fmean-fmin),100*(fmax-fmean)); 
%er = errorbar(1:8,100*fmean,100*(fmin),100*(fmax)); 
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca,'XTickLabel',names)
ylabel('Percent change')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_all_types_ensem_minmax.png'])

%% size
sbar = pbar(1:2);
sbar(3:4) = fmean(1:2);
smean = fmean(1:2);
sstd = fstd(1:2);
smax = fmax(1:2);
smin = fmin(1:2);

figure(5)
bar(sbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),100*sstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','Mesozoo','M Fishes','L Fishes'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_ensem_1std_color.png'])

figure(6)
bar(sbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),200*sstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','Mesozoo','M Fishes','L Fishes'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_ensem_2std.png'])

figure(7)
bar(sbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),100*(smean-smin),100*(smax-smean),'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','Mesozoo','M Fishes','L Fishes'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_ensem_minmax.png'])

%% BW
figure(15)
b=bar(sbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),100*sstd,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','Mesozoo','M Fishes','L Fishes'})
ylabel('Percent change')
%title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_ensem_1std.png'])

%% type
tbar = fmean(3:5);
tstd = fstd(3:5);
tmax = fmax(3:5);
tmin = fmin(3:5);

figure(8)
bar(tbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(1:3,100*tbar,100*tstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_type_ensem_1std.png'])

figure(9)
bar(tbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(1:3,100*tbar,200*tstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_type_ensem_2std.png'])

figure(10)
bar(tbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 0])
er = errorbar(1:3,100*tbar,100*(tbar-tmin),100*(tmax-tbar),'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_type_ensem_minmax.png'])

%% type with bent
bbar = fmean([3:5 8]);
bstd = fstd([3:5 8]);

figure(11)
bar(bbar*100,'k'); hold on;
%colormap('gray')
ylim([-25 10])
er = errorbar(1:4,100*bbar,100*bstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D','B'})
ylabel('Percent change')
title('Global change in productivity')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_typeB_ensem_1std.png'])

%% Type & size subplot
figure(12)
subplot(2,2,1)
bar(sbar*100,'k'); hold on;
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),100*sstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','MesoZ','M','L'})
ylabel('Percent change')

subplot(2,2,2)
bar(tbar*100,'k'); hold on;
ylim([-25 0])
er = errorbar(1:3,100*tbar,100*tstd,'LineWidth',2);    
er.Color = [0 0.5 1];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D'})
%ylabel('Percent change')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_type_ensem_1std.png'])

%%
figure(13)
subplot(2,2,1)
b=bar(pbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 10])
set(gca,'XTickLabel',{'NPP','MesoZ','Det'})
ylabel('Percent change')

subplot(2,2,2)
b=bar(bbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 10])
er = errorbar(1:4,100*bbar,100*bstd,'LineWidth',1.5);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D','B'})
%ylabel('Percent change')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_bgc_type_ensem_1std.png'])

%% tables
tstats(:,1) = fmean';
tstats(:,2) = fstd';
tstats(:,3) = fmin';
tstats(:,4) = fmax';
Stab = array2table(tstats,'VariableNames',{'mean','std','min','max'},...
    'RowNames',names);
writetable(Stab,[nfile 'Hist_Fore_All_fish03_ensem6_mid_kt3_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)







