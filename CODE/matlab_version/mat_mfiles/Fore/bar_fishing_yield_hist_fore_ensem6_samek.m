% Bar plot of changes in fishing yield
% with std dev from ensemble sims
% 5 most sens params plus same k's

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
efile = 'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',efile,'/full_runs/'];
dfile = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];

%% Original parameters
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
fpath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

%% Historic & Forecast
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
    'HF_tsyc','HP_tsyc','HD_tsyc','HA_tsyc','HMF_tsyc','HLP_tsyc',...
    'HLD_tsyc','HAA_tsyc','FF_tsyc','FP_tsyc','FD_tsyc','FA_tsyc',...
    'FMF_tsyc','FLP_tsyc','FLD_tsyc','FAA_tsyc');

%% Ensemble parameter sets
epath = nfile;
epath2 = dfile;

load([epath 'Historic_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'Hlme_tsc_mf','Hlme_tsc_mp','Hlme_tsc_md','Hlme_tsc_lp','Hlme_tsc_ld');
load([epath 'Forecast_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'Flme_tsc_mf','Flme_tsc_mp','Flme_tsc_md','Flme_tsc_lp','Flme_tsc_ld');

% Time series of catch (g per km2 per year)

%% ts
%In original saved file
y1 = 1861:2005;
y2 = 2006:2100;
y = [y1 y2];

HF = Hlme_tsc_mf;
HP = Hlme_tsc_mp + Hlme_tsc_lp;
HD = Hlme_tsc_md + Hlme_tsc_ld;
HLP = Hlme_tsc_lp;
HLD = Hlme_tsc_ld;

%%
HA = HF + HP + HD;
HAA = HF + HLP + HLD;

FF = Flme_tsc_mf;
FP = Flme_tsc_mp + Flme_tsc_lp;
FD = Flme_tsc_md + Flme_tsc_ld;
FA = FF + FP + FD;
FLP = Flme_tsc_lp;
FLD = Flme_tsc_ld;
FAA = FF + FLP + FLD;

tForig = [HF_tsyc FF_tsyc];
tPorig = [HP_tsyc FP_tsyc];
tDorig = [HD_tsyc FD_tsyc];
tAorig = [HA_tsyc FA_tsyc];
tLPorig = [HLP_tsyc FLP_tsyc];
tLDorig = [HLD_tsyc FLD_tsyc];
tAAorig = [HAA_tsyc FAA_tsyc];

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];
tLP = [HLP FLP];
tLD = [HLD FLD];
tAA = [HAA FAA];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];
tLP = [tLP; tLPorig];
tLD = [tLD; tLDorig];
tAA = [tAA; tAAorig];

%by trophic level to compare to Moore et al 2018
t3 = tF;
t4 = tLP + tLD;

%% last 50 yrs of each
h50 = find(y>1950 & y<=2000);
f50 = find(y>2050 & y<=2100);

fbar(:,1) = ((mean(tF(:,f50),2))-(mean(tF(:,h50),2))) ./ (mean(tF(:,h50),2));
fbar(:,2) = ((mean(tP(:,f50),2))-(mean(tP(:,h50),2))) ./ (mean(tP(:,h50),2));
fbar(:,3) = ((mean(tD(:,f50),2))-(mean(tD(:,h50),2))) ./ (mean(tD(:,h50),2));
fbar(:,4) = ((mean(tA(:,f50),2))-(mean(tA(:,h50),2))) ./ (mean(tA(:,h50),2));
fbar(:,5) = ((mean(tLP(:,f50),2))-(mean(tLP(:,h50),2))) ./ (mean(tLP(:,h50),2));
fbar(:,6) = ((mean(tLD(:,f50),2))-(mean(tLD(:,h50),2))) ./ (mean(tLD(:,h50),2));
fbar(:,7) = ((mean(tAA(:,f50),2))-(mean(tAA(:,h50),2))) ./ (mean(tAA(:,h50),2));
fbar(:,8) = ((mean(t3(:,f50),2))-(mean(t3(:,h50),2))) ./ (mean(t3(:,h50),2));
fbar(:,9) = ((mean(t4(:,f50),2))-(mean(t4(:,h50),2))) ./ (mean(t4(:,h50),2));
names = {'F','P','D','All','LP','LD','Adults','TL3','TL4'};

dbar(:,1) = ((mean(tF(:,f50),2))-(mean(tF(:,h50),2))) ;
dbar(:,2) = ((mean(tP(:,f50),2))-(mean(tP(:,h50),2))) ;
dbar(:,3) = ((mean(tD(:,f50),2))-(mean(tD(:,h50),2))) ;
dbar(:,4) = ((mean(tA(:,f50),2))-(mean(tA(:,h50),2))) ;
dbar(:,5) = ((mean(tLP(:,f50),2))-(mean(tLP(:,h50),2))) ;
dbar(:,6) = ((mean(tLD(:,f50),2))-(mean(tLD(:,h50),2))) ;
dbar(:,7) = ((mean(tAA(:,f50),2))-(mean(tAA(:,h50),2))) ;
dbar(:,8) = ((mean(t3(:,f50),2))-(mean(t3(:,h50),2))) ;
dbar(:,9) = ((mean(t4(:,f50),2))-(mean(t4(:,h50),2))) ;

%% take mean and error bars
%means of global means, should they be means of all cells???
fmean = mean(fbar);
fstd = std(fbar);
fmax = max(fbar);
fmin = min(fbar);

abar  = fbar([1,5:7]);
amean = fmean([1,5:7]);
astd  = fstd([1,5:7]);

%% bar graphs
% 1 st dev
figure(1)
b=bar(100*fmean); hold on;
b.FaceColor = [0 0.5 0.75];
er = errorbar(1:9,100*fmean,100*fstd);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca,'XTickLabel',names)
ylabel('Percent change')
print('-dpng',[ppath 'Hist_Fore_All_fish03_bar_yield_all_types_ensem_samek_1std.png'])

%% BW
% All
figure(2)
b=bar(fmean(1:4)*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(1:4,100*fmean(1:4),100*fstd(1:4),'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',names(1:4))
ylabel('Percent change')
%title('Global change in productivity')
print('-dpng',[ppath 'Hist_Fore_All_fish03_bar_yield_types_ensem_samek_1std.png'])

%% Adults
figure(3)
b=bar(amean*100,'k'); hold on;
b.FaceColor = [0.65 0.65 0.65];
ylim([-25 0])
er = errorbar(1:4,100*amean,100*astd,'LineWidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'XTickLabel',{'MF','LP','LD','All adults'})
ylabel('Percent change')
%title('Global change in productivity')
print('-dpng',[ppath 'Hist_Fore_All_fish03_bar_yield_adults_ensem_samek_1std.png'])

%% tables
tstats(:,1) = fmean';
tstats(:,2) = fstd';
tstats(:,3) = fmin';
tstats(:,4) = fmax';
Stab = array2table(tstats,'VariableNames',{'mean','std','min','max'},...
    'RowNames',names);
writetable(Stab,[epath 'Hist_Fore_All_fish03_ensem6_mid_samek_yield_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Stab,[epath2 'Hist_Fore_All_fish03_ensem6_mid_samek_yield_pdiff_MeanStd.csv'],...
    'Delimiter',',','WriteRowNames',true)

save([epath 'Hist_Fore_All_fish03_ensem6_mid_samek_yield_diff_pdiff.mat'],...
    'fbar','dbar','names')
save([epath2 'Hist_Fore_All_fish03_ensem6_mid_samek_yield_diff_pdiff.mat'],...
    'fbar','dbar','names')







