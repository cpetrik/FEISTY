% Plot Change in Prod vs Nu
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Ensemble mid6, temp3
% Ensemble end members with max and min change in production relative 
% to their 1951 value
% Nu scaled by Cmax

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');

%% Temp & Det
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'ptemp_mean_hist',...
    'ptemp_mean_fore','btemp_mean_hist','btemp_mean_fore','det_mean_hist','det_mean_fore');

grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

ptemp_hist = ptemp_mean_hist(ID);
ptemp_fore = ptemp_mean_fore(ID);
btemp_hist = btemp_mean_hist(ID);
btemp_fore = btemp_mean_fore(ID);

det_hist = det_mean_hist(ID) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore(ID) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

%% FEISTY Output
% cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
% fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

%% Ensemble parameter sets
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

% Production
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'hSPmF','hSPmP','hSPmD','hSPlP','hSPlD','hSPB')
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'fSPmF','fSPmP','fSPmD','fSPlP','fSPlD','fSPB')

% MaxMin Prod sims
load([epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat'],...
    'Ttab','tstats','Ptab','pstats');

% Nu 
load([epath 'Historic_All_fish03_ensem6_MaxMin1951_nus.mat'],...
    'hNSmF','hNSmP','hNSmD','hNSlP','hNSlD')
load([epath 'Forecast_All_fish03_ensem6_MaxMin1951_nus.mat'],...
    'fNSmF','fNSmP','fNSmD','fNSlP','fNSlD',...
    'fnms','snms');

load([epath 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext');
params = red_params;
psets = unique(pstats(:,1));
nparam = length(psets);

%% calculate: 
%benthic production = benthic biomass * detritus flux * benthic efficiency
hPB = hSPB .* det_hist .* 0.075;
fPB = fSPB .* det_fore .* 0.075;

%% Calc cmax
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

HMPcmax = (exp(0.063.*(ptemp_hist-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
HLPcmax = (exp(0.063.*(ptemp_hist-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;
HMBcmax = (exp(0.063.*(btemp_hist-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
HLBcmax = (exp(0.063.*(btemp_hist-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;

FMPcmax = (exp(0.063.*(ptemp_fore-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
FLPcmax = (exp(0.063.*(ptemp_fore-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;
FMBcmax = (exp(0.063.*(btemp_fore-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
FLBcmax = (exp(0.063.*(btemp_fore-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;

% 
HsnuMF = hNSmF ./ repmat(HMPcmax,1,nparam);
HsnuMP = hNSmP ./ repmat(HMPcmax,1,nparam);
HsnuMD = hNSmD ./ repmat(HMBcmax,1,nparam);
HsnuLP = hNSlP ./ repmat(HLPcmax,1,nparam);
HsnuLD = hNSlD ./ repmat(HLBcmax,1,nparam);

FsnuMF = fNSmF ./ repmat(FMPcmax,1,nparam);
FsnuMP = fNSmP ./ repmat(FMPcmax,1,nparam);
FsnuMD = fNSmD ./ repmat(FMBcmax,1,nparam);
FsnuLP = fNSlP ./ repmat(FLPcmax,1,nparam);
FsnuLD = fNSlD ./ repmat(FLBcmax,1,nparam);

%% Reduce prod to just MaxMinDiff sims
hPmF = hSPmF(:,psets);
hPmP = hSPmP(:,psets);
hPmD = hSPmD(:,psets);
hPlP = hSPlP(:,psets);
hPlD = hSPlD(:,psets);
hPB  = hPB(:,psets);

fPmF = fSPmF(:,psets);
fPmP = fSPmP(:,psets);
fPmD = fSPmD(:,psets);
fPlP = fSPlP(:,psets);
fPlD = fSPlD(:,psets);
fPB  = fPB(:,psets);

%% Calc differences
dPmF = (fPmF - hPmF);
dPmP = (fPmP - hPmP);
dPmD = (fPmD - hPmD);
dPlP = (fPlP - hPlP);
dPlD = (fPlD - hPlD);
dPB  = (fPB  - hPB);

dNmF = (FsnuMF - HsnuMF);
dNmP = (FsnuMP - HsnuMP);
dNmD = (FsnuMD - HsnuMD);
dNlP = (FsnuLP - HsnuLP);
dNlD = (FsnuLD - HsnuLD);

%% MF max
xmf = find(psets==pstats(3,1));
figure(1)
subplot(2,2,1)
scatter(HsnuMF(:,xmf), hPmF(:,xmf), 10, ptemp_hist, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod MF')
title('Hist')
axis([-0.05 0.45 0 0.3])

subplot(2,2,2)
scatter(FsnuMF(:,xmf), fPmF(:,xmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod MF')
title('Fore')
axis([-0.05 0.45 0 0.3])

subplot(2,2,3)
scatter(FsnuMF(:,xmf), dPmF(:,xmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax Fore')
ylabel('\Delta Prod MF')
title('max \Delta Prod F')
axis([-0.05 0.45 -0.065 0.06])

subplot(2,2,4)
scatter(dNmF(:,xmf), dPmF(:,xmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod MF')
title('max \Delta Prod F')
axis([-0.15 0.25 -0.065 0.06])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_scnu_vs_prod_maxF.png'])

%% MF min
nmf = find(psets==pstats(4,1));
figure(2)
subplot(2,2,1)
scatter(HsnuMF(:,nmf), hPmF(:,nmf), 10, ptemp_hist, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod MF')
title('Hist')
axis([-0.05 0.45 0 0.3])

subplot(2,2,2)
scatter(FsnuMF(:,nmf), fPmF(:,nmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod MF')
title('Fore')
axis([-0.05 0.45 0 0.3])

subplot(2,2,3)
scatter(FsnuMF(:,nmf), dPmF(:,nmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax Fore')
ylabel('\Delta Prod MF')
title('min \Delta Prod F')
axis([-0.05 0.45 -0.06 0.06])

subplot(2,2,4)
scatter(dNmF(:,nmf), dPmF(:,nmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod MF')
title('min \Delta Prod F')
axis([-0.15 0.25 -0.065 0.06])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_scnu_vs_prod_minF.png'])


%% LP max
xp = find(psets==pstats(5,1));
figure(3)
subplot(2,2,1)
scatter(HsnuLP(:,xp), hPlP(:,xp), 10, ptemp_hist, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LP')
title('Hist')
axis([-0.5 0.3 0 0.05])

subplot(2,2,2)
scatter(FsnuLP(:,xp), fPlP(:,xp), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LP')
title('Fore')
axis([-0.5 0.3 0 0.05])

subplot(2,2,3)
scatter(FsnuLP(:,xp), dPlP(:,xp), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax Fore')
ylabel('\Delta Prod P')
title('max \Delta Prod P')
axis([-0.5 0.3 -0.05 0.02])

subplot(2,2,4)
scatter(dNlP(:,xp), dPlP(:,xp), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod P')
title('max \Delta Prod P')
axis([-0.3 0.5 -0.05 0.02])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_scnu_vs_prod_maxP.png'])

%% LP min
np = find(psets==pstats(6,1));
figure(4)
subplot(2,2,1)
scatter(HsnuLP(:,np), hPlP(:,np), 10, ptemp_hist, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LP')
title('Hist')
axis([-0.5 0.3 0 0.05])

subplot(2,2,2)
scatter(FsnuLP(:,np), fPlP(:,np), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LP')
title('Fore')
axis([-0.5 0.3 0 0.05])

subplot(2,2,3)
scatter(FsnuLP(:,np), dPlP(:,np), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax Fore')
ylabel('\Delta Prod P')
title('min \Delta Prod P')
axis([-0.5 0.3 -0.05 0.02])

subplot(2,2,4)
scatter(dNlP(:,np), dPlP(:,np), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod P')
title('min \Delta Prod P')
axis([-0.3 0.5 -0.05 0.02])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_scnu_vs_prod_minP.png'])


%% LD max
xd = find(psets==pstats(7,1));
figure(5)
subplot(2,2,1)
scatter(HsnuLD(:,xd), hPlD(:,xd), 10, btemp_hist, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LD')
title('Hist')
axis([-0.05 0.3 0 0.05])

subplot(2,2,2)
scatter(FsnuLD(:,xd), fPlD(:,xd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LD')
title('Fore')
axis([-0.05 0.3 0 0.05])

subplot(2,2,3)
scatter(FsnuLD(:,xd), dPlD(:,xd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax Fore')
ylabel('\Delta Prod D')
title('max \Delta Prod D')
axis([-0.05 0.3 -0.015 0.01])

subplot(2,2,4)
scatter(dNlD(:,xd), dPlD(:,xd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod D')
title('max \Delta Prod D')
axis([-0.07 0.08 -0.015 0.01])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_scnu_vs_prod_maxD.png'])

%% LD min
nd = find(psets==pstats(8,1));
figure(6)
subplot(2,2,1)
scatter(HsnuLD(:,nd), hPlD(:,nd), 10, btemp_hist, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LD')
title('Hist')
axis([-0.05 0.3 0 0.05])

subplot(2,2,2)
scatter(FsnuLD(:,nd), fPlD(:,nd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax')
ylabel('Prod LD')
title('Fore')
axis([-0.05 0.3 0 0.05])

subplot(2,2,3)
scatter(FsnuLD(:,nd), dPlD(:,nd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('Nu/Cmax Fore')
ylabel('\Delta Prod D')
title('min \Delta Prod D')
axis([-0.05 0.3 -0.015 0.01])

subplot(2,2,4)
scatter(dNlD(:,nd), dPlD(:,nd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod D')
title('min \Delta Prod D')
axis([-0.07 0.08 -0.015 0.01])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_scnu_vs_prod_minD.png'])

%%
subplot(3,2,1)
scatter(ptemp_hist, hPmF(:,nd), 10,  'filled')
% cmocean('thermal')
% colorbar
% caxis([0 30])
xlabel('temp')
ylabel('MF prod')
title('Hist')
%axis([-0.05 0.3 -0.015 0.01])

subplot(3,2,2)
scatter(ptemp_fore, fPmF(:,nd), 10,  'filled')
% cmocean('thermal')
% colorbar
% caxis([0 30])
xlabel('temp')
ylabel('MF prod')
title('Fore')

subplot(3,2,3)
scatter(ptemp_hist, hPlP(:,nd), 10,  'filled')
% cmocean('thermal')
% colorbar
% caxis([0 30])
xlabel('temp')
ylabel('LP prod')
title('Hist')
%axis([-0.05 0.3 -0.015 0.01])

subplot(3,2,4)
scatter(ptemp_fore, fPlP(:,nd), 10,  'filled')
% cmocean('thermal')
% colorbar
% caxis([0 30])
xlabel('temp')
ylabel('LP prod')
title('Fore')

subplot(3,2,5)
scatter(btemp_hist, hPlD(:,nd), 10,  'filled')
% cmocean('thermal')
% colorbar
% caxis([0 30])
xlabel('temp')
ylabel('LD prod')
title('Hist')
%axis([-0.05 0.3 -0.015 0.01])

subplot(3,2,6)
scatter(btemp_fore, fPlD(:,nd), 10,  'filled')
% cmocean('thermal')
% colorbar
% caxis([0 30])
xlabel('temp')
ylabel('LD prod')
title('Fore')

