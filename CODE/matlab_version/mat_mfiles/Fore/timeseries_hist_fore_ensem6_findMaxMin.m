% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Original parameters
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];
load([fpath 'Time_Means_Historic_Forecast_',harv,'_' cfile '.mat']);

%% Ensemble parameter sets
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat']);

lpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];

%% ts
%In original saved file
% y1 = 1860+(1/12):(1/12):2005;
% y2 = 2005+(1/12):(1/12):2100;
% y = [y1 y2];

% SOMETHING WRONG WITH MOVING MEAN OF #28 AROUND 1949
% FIX TS OF SMALL FISH AT T=1080
%tF[1081,28]=1.94e+35
%tP[1081,28]=1.94e+35
%tD[1081,28]=1.94e+35
%tA[1081,28]=5.83e+35
hTsF(28,1081) = (hTsF(28,1080) + hTsF(28,1082))/2;
hTsP(28,1081) = (hTsP(28,1080) + hTsP(28,1082))/2;
hTsD(28,1081) = (hTsD(28,1080) + hTsD(28,1082))/2;

HF = hTsF + hTmF;
HP = hTsP + hTmP + hTlP;
HD = hTsD + hTmD + hTlD;
HA = HF + HP + HD;

FF = fTsF + fTmF;
FP = fTsP + fTmP + fTlP;
FD = fTsD + fTmD + fTlD;
FA = FF + FP + FD;

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];

%% Calc variability at 1900, 2000, and 2100

%Table with most & least fish
tabML(1,1) = find(FA(:,1140)==max(FA(:,1140)));
tabML(1,2) = find(FA(:,1140)==min(FA(:,1140)));
tabML(2,1) = find(FF(:,1140)==max(FF(:,1140)));
tabML(2,2) = find(FF(:,1140)==min(FF(:,1140)));
tabML(3,1) = find(FP(:,1140)==max(FP(:,1140)));
tabML(3,2) = find(FP(:,1140)==min(FP(:,1140)));
tabML(4,1) = find(FD(:,1140)==max(FD(:,1140)));
tabML(4,2) = find(FD(:,1140)==min(FD(:,1140)));

Ftab = array2table(tabML,'VariableNames',{'Most','Least'},...
    'RowNames',{'All Fish','F','P','D'});
writetable(Ftab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MostLeast.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Ftab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_MostLeast.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% Variability and difference
tstats(1,1) = var(tA(:,480));
tstats(2,1) = var(tA(:,1680));
tstats(3,1) = var(tA(:,2880));

tstats(1,2) = max(tA(:,480)) - min(tA(:,480));
tstats(2,2) = max(tA(:,1680)) - min(tA(:,1680));
tstats(3,2) = max(tA(:,2880)) - min(tA(:,2880));

tstats(1,3) = var(tF(:,480));
tstats(2,3) = var(tF(:,1680));
tstats(3,3) = var(tF(:,2880));

tstats(1,4) = max(tF(:,480)) - min(tF(:,480));
tstats(2,4) = max(tF(:,1680)) - min(tF(:,1680));
tstats(3,4) = max(tF(:,2880)) - min(tF(:,2880));

tstats(1,5) = var(tP(:,480));
tstats(2,5) = var(tP(:,1680));
tstats(3,5) = var(tP(:,2880));

tstats(1,6) = max(tP(:,480)) - min(tP(:,480));
tstats(2,6) = max(tP(:,1680)) - min(tP(:,1680));
tstats(3,6) = max(tP(:,2880)) - min(tP(:,2880));

tstats(1,7) = var(tD(:,480));
tstats(2,7) = var(tD(:,1680));
tstats(3,7) = var(tD(:,2880));

tstats(1,8) = max(tD(:,480)) - min(tD(:,480));
tstats(2,8) = max(tD(:,1680)) - min(tD(:,1680));
tstats(3,8) = max(tD(:,2880)) - min(tD(:,2880));

Stab = array2table(tstats,'VariableNames',{'varAll','diffAll','varF','diffF',...
    'varP','diffP','varD','diffD'},...
    'RowNames',{'1900','2000','2100'});
writetable(Stab,[epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarDiff.csv'],...
    'Delimiter',',','WriteRowNames',true)
writetable(Stab,[lpath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarDiff.csv'],...
    'Delimiter',',','WriteRowNames',true)

%% Line color order
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

%% moving means
mmtF = movmean(tF,31,2); %[1081,28]=1.94e+35
mmtP = movmean(tP,31,2); %[1081,28]=1.94e+35
mmtD = movmean(tD,31,2); %[1081,28]=1.94e+35
mmtA = movmean(tA,31,2); %[1081,28]=5.83e+35

mmoF = movmean(tForig,31);
mmoP = movmean(tPorig,31);
mmoD = movmean(tDorig,31);
mmoA = movmean(tAorig,31);

%% Individual Plots 
% All
figure(1)
subplot(3,3,1)
plot(y,log10(mmtA(1:5,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('1','2','3','4','5')
legend('location','southwest')

subplot(3,3,2)
plot(y,log10(mmtA(6:10,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('6','7','8','9','10')
legend('location','southwest')
title('All')

subplot(3,3,3)
plot(y,log10(mmtA(11:15,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('11','12','13','14','15')
legend('location','southwest')

subplot(3,3,4)
plot(y,log10(mmtA(16:20,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('16','17','18','19','20')
legend('location','southwest')

subplot(3,3,5)
plot(y,log10(mmtA(21:25,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('21','22','23','24','25')
legend('location','southwest')

subplot(3,3,6)
plot(y,log10(mmtA(26:30,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('26','27','28','29','30')
legend('location','southwest')

subplot(3,3,7)
plot(y,log10(mmtA(31:35,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('31','32','33','34','35')
legend('location','southwest')

subplot(3,3,8)
plot(y,log10(mmtA(36:40,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('36','37','38','39','40')
legend('location','southwest')

subplot(3,3,9)
plot(y,log10(mmtA(41:43,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend('41','42','43')
legend('location','southwest')

print('-dpng',[ppath 'Hist_Fore_',harv,'_all_types_ensem_mid6_temp3_lab.png'])

%% F
figure(2)
subplot(3,3,1)
plot(y,log10(mmtF(1:5,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('1','2','3','4','5')
legend('location','southwest')

subplot(3,3,2)
plot(y,log10(mmtF(6:10,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('6','7','8','9','10')
legend('location','southwest')
title('Forage')

subplot(3,3,3)
plot(y,log10(mmtF(11:15,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('11','12','13','14','15')
legend('location','southwest')

subplot(3,3,4)
plot(y,log10(mmtF(16:20,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('16','17','18','19','20')
legend('location','southwest')

subplot(3,3,5)
plot(y,log10(mmtF(21:25,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('21','22','23','24','25')
legend('location','southwest')

subplot(3,3,6)
plot(y,log10(mmtF(26:30,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('26','27','28','29','30')
legend('location','southwest')

subplot(3,3,7)
plot(y,log10(mmtF(31:35,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('31','32','33','34','35')
legend('location','southwest')

subplot(3,3,8)
plot(y,log10(mmtF(36:40,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('36','37','38','39','40')
legend('location','southwest')

subplot(3,3,9)
plot(y,log10(mmtF(41:43,:))); hold on;
xlim([1900 2100])
ylim([0.05 0.35])
legend('41','42','43')
legend('location','southwest')
print('-dpng',[ppath 'Hist_Fore_',harv,'_forage_ensem_mid6_temp3_lab.png'])

%% P
figure(3)
subplot(3,3,1)
plot(y,log10(mmtP(1:5,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('1','2','3','4','5')
legend('location','southwest')

subplot(3,3,2)
plot(y,log10(mmtP(6:10,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('6','7','8','9','10')
legend('location','southwest')
title('Large pelagics')

subplot(3,3,3)
plot(y,log10(mmtP(11:15,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('11','12','13','14','15')
legend('location','southwest')

subplot(3,3,4)
plot(y,log10(mmtP(16:20,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('16','17','18','19','20')
legend('location','southwest')

subplot(3,3,5)
plot(y,log10(mmtP(21:25,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('21','22','23','24','25')
legend('location','southwest')

subplot(3,3,6)
plot(y,log10(mmtP(26:30,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('26','27','28','29','30')
legend('location','southwest')

subplot(3,3,7)
plot(y,log10(mmtP(31:35,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('31','32','33','34','35')
legend('location','southwest')

subplot(3,3,8)
plot(y,log10(mmtP(36:40,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('36','37','38','39','40')
legend('location','southwest')

subplot(3,3,9)
plot(y,log10(mmtP(41:43,:))); hold on;
xlim([1900 2100])
ylim([-0.2 0.4])
legend('41','42','43')
legend('location','southwest')
print('-dpng',[ppath 'Hist_Fore_',harv,'_pel_ensem_mid6_temp3_lab.png'])

%% D
figure(4)
subplot(3,3,1)
plot(y,log10(mmtD(1:5,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('1','2','3','4','5')
legend('location','southwest')

subplot(3,3,2)
plot(y,log10(mmtD(6:10,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('6','7','8','9','10')
legend('location','southwest')
title('Demersals')

subplot(3,3,3)
plot(y,log10(mmtD(11:15,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('11','12','13','14','15')
legend('location','southwest')

subplot(3,3,4)
plot(y,log10(mmtD(16:20,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('16','17','18','19','20')
legend('location','southwest')

subplot(3,3,5)
plot(y,log10(mmtD(21:25,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('21','22','23','24','25')
legend('location','southwest')

subplot(3,3,6)
plot(y,log10(mmtD(26:30,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('26','27','28','29','30')
legend('location','southwest')

subplot(3,3,7)
plot(y,log10(mmtD(31:35,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('31','32','33','34','35')
legend('location','southwest')

subplot(3,3,8)
plot(y,log10(mmtD(36:40,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('36','37','38','39','40')
legend('location','southwest')

subplot(3,3,9)
plot(y,log10(mmtD(41:43,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend('41','42','43')
legend('location','southwest')
print('-dpng',[ppath 'Hist_Fore_',harv,'_dem_ensem_mid6_temp3_lab.png'])


%%
% Line color order
cm5=[0 0.5 0;... %dk green
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 0;...     %r
    0 0 0]; 

set(groot,'defaultAxesColorOrder',cm5);

spec = [15;17;28;32;39];
ltext = {'min \DeltaD','min \DeltaMLFPAll','max \DeltaMFAll','max \DeltaD','max \DeltaP'};

figure(5)
subplot(2,2,4)
plot(y,log10(mmtA(spec,:))); hold on;
xlim([1900 2100])
ylim([0.425 0.725])
legend(ltext)
legend('location','southwest')
title('All')

subplot(2,2,1)
plot(y,log10(mmtF(spec,:))); hold on;
xlim([1900 2100])
ylim([0.0 0.35])
% legend(num2str(spec))
% legend('location','southwest')
title('Forage')

subplot(2,2,2)
plot(y,log10(mmtP(spec,:))); hold on;
xlim([1900 2100])
ylim([-0.1 0.3])
legend(ltext)
legend('location','southwest')
title('Large pelagics')

subplot(2,2,3)
plot(y,log10(mmtD(spec,:))); hold on;
xlim([1900 2100])
ylim([-0.15 0.1])
legend(ltext)
legend('location','southwest')
title('Demersals')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_lab_maxmin_ts.png'])



