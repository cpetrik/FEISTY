% Bar plot of trophic amplification
% with std dev from ensemble sims var temp and equal temp
% NO benthic production = benthic biomass * detritus flux * benthic efficiency
% calc benthic production = detritus flux * benthic efficiency

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

%% Varying temp-dep
% epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
%     'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
% load([epath 'Prod_diff_50yr_ensem6_mid_kt2_bestAIC_multFup_multPneg.mat'])
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([nfile 'Prod_diff_50yr_ensem6_mid_kt3_bestAIC_multFup_multPneg.mat'])

% take mean and error bars
%functional groups and sizes
fmean = mean(fbar);
fstd = std(fbar);
fmax = max(fbar);
fmin = min(fbar);

%size
sbar = pbar(1:2);
sbar(3:4) = fmean(1:2);
smean = fmean(1:2);
sstd = fstd(1:2);
smax = fmax(1:2);
smin = fmin(1:2);

%type
tbar = fmean(3:5);
tstd = fstd(3:5);
tmax = fmax(3:5);
tmin = fmin(3:5);

%type with bent
bbar = fmean([3:5 8]);
bstd = fstd([3:5 8]);

clear fbar 

%% Equal temp dep
efile = 'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';
dfile = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
load([nfile 'Prod_diff_50yr_ensem6_mid_samek_bestAIC_Fupneg_mult8_Pneg2_mult3.mat'])

% take mean and error bars
fmean2 = mean(fbar);
fstd2 = std(fbar);
fmax2 = max(fbar);
fmin2 = min(fbar);

sbar2 = pbar(1:2);
sbar2(3:4) = fmean2(1:2);
smean2 = fmean2(1:2);
sstd2 = fstd2(1:2);

tbar2 = fmean2(3:5);
tstd2 = fstd2(3:5);
tmax2 = fmax2(3:5);
tmin2 = fmin2(3:5);

bbar2 = fmean2([3:5 8]);
bstd2 = fstd2([3:5 8]);

%% bar graphs
% BW
figure(1)
subplot(2,2,1)
b=bar(sbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(3:4,100*sbar(3:4),100*sstd,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','MesoZ','M Fish','L Fish'})
ylabel('Percent change')
text(0,1,'A')

subplot(2,2,2)
b=bar(sbar2*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(3:4,100*sbar2(3:4),100*sstd2,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'NPP','MesoZ','M Fish','L Fish'})
text(0,1,'B')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_size_both_ensem_1std.png'])

%% Type & size subplot
figure(2)
subplot(3,3,1)
b=bar(pbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
set(gca,'XTickLabel',{'NPP','MesoZ','Det'})
ylabel('Percent change')
text(0,2,'A')

subplot(3,3,2)
b=bar(bbar*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(1:3,100*bbar(1:3),100*bstd(1:3),'LineWidth',1.5);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D','B'})
text(0,2,'B')

subplot(3,3,3)
b=bar(bbar2*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(1:3,100*bbar2(1:3),100*bstd2(1:3),'LineWidth',1.5);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'F','P','D','B'})
text(0,2,'C')
print('-dpng',[pp 'Hist_Fore_All_fish03_troph_amp_bgc_type_both_ensem_1std.png'])








