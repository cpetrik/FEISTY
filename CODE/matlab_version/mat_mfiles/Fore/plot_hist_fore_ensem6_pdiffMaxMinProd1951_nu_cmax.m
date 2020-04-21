% Plot Change in Prod vs Nu
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Ensemble mid6, temp3
% Ensemble end members with max and min percent change in production relative 
% to their 1951 value
% Nu scaled by Cmax (both FEISTY output)

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
    'Mtab','mmstat','Dtab','pdstat');

% Nu 
load([epath 'Historic_All_fish03_ensem6_pdiffMaxMin1951_nus.mat'],...
    'hNSmF','hNSmP','hNSmD','hNSlP','hNSlD')
load([epath 'Forecast_All_fish03_ensem6_pdiffMaxMin1951_nus.mat'],...
    'fNSmF','fNSmP','fNSmD','fNSlP','fNSlD',...
    'fnms','snms');

% Cmax
load([epath 'Historic_All_fish03_ensem6_pdiffMaxMin1951_cmaxs.mat'],...
    'hCSsF','hCSsP','hCSsD','hCSmF','hCSmP','hCSmD','hCSlP','hCSlD')
load([epath 'Forecast_All_fish03_ensem6_pdiffMaxMin1951_cmaxs.mat'],...
    'fCSsF','fCSsP','fCSsD','fCSmF','fCSmP','fCSmD','fCSlP','fCSlD');

load([epath 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext');
params = red_params;
psets = unique(pdstat(:,1));
nparam = length(psets);

%% calculate: 
%benthic production = benthic biomass * detritus flux * benthic efficiency
hPB = hSPB .* det_hist .* 0.075;
fPB = fSPB .* det_fore .* 0.075;

%% Scale by cmax 
HsnuMF = hNSmF ./ hCSmF;
HsnuMP = hNSmP ./ hCSmP;
HsnuMD = hNSmD ./ hCSmD;
HsnuLP = hNSlP ./ hCSlP;
HsnuLD = hNSlD ./ hCSlD;

FsnuMF = fNSmF ./ fCSmF;
FsnuMP = fNSmP ./ fCSmP;
FsnuMD = fNSmD ./ fCSmD;
FsnuLP = fNSlP ./ fCSlP;
FsnuLD = fNSlD ./ fCSlD;

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
dPmF = (fPmF - hPmF) ;
dPmP = (fPmP - hPmP) ;
dPmD = (fPmD - hPmD) ;
dPlP = (fPlP - hPlP) ;
dPlD = (fPlD - hPlD) ;
dPB  = (fPB  - hPB) ;

pdPmF = (fPmF - hPmF) ./hPmF;
pdPmP = (fPmP - hPmP) ./hPmP;
pdPmD = (fPmD - hPmD) ./hPmD;
pdPlP = (fPlP - hPlP) ./hPlP;
pdPlD = (fPlD - hPlD) ./hPlD;
pdPB  = (fPB  - hPB) ./hPB;

dNmF = (FsnuMF - HsnuMF);
dNmP = (FsnuMP - HsnuMP);
dNmD = (FsnuMD - HsnuMD);
dNlP = (FsnuLP - HsnuLP);
dNlD = (FsnuLD - HsnuLD);

%% MF max
xmf = find(psets==pdstat(3,1));
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
title('max % \Delta Prod F')
axis([-0.05 0.45 -0.065 0.06])

subplot(2,2,4)
scatter(dNmF(:,xmf), dPmF(:,xmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod MF')
title('max % \Delta Prod F')
axis([-0.15 0.25 -0.065 0.06])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_nucmax_vs_prod_maxF_pdiff.png'])

%% MF min
nmf = find(psets==pdstat(4,1));
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
title('min % \Delta Prod F')
axis([-0.05 0.45 -0.065 0.06])

subplot(2,2,4)
scatter(dNmF(:,nmf), dPmF(:,nmf), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('\Delta Prod MF')
title('min % \Delta Prod F')
axis([-0.15 0.25 -0.065 0.06])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_nucmax_vs_prod_minF_pdiff.png'])


%% LP max
xp = find(psets==pdstat(5,1));
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
ylabel('% \Delta Prod P')
title('max % \Delta Prod P')
axis([-0.5 0.3 -0.05 0.02])

subplot(2,2,4)
scatter(dNlP(:,xp), dPlP(:,xp), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('% \Delta Prod P')
title('max % \Delta Prod P')
axis([-0.3 0.5 -0.05 0.02])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_nucmax_vs_prod_maxP_pdiff.png'])

%% LP min
np = find(psets==pdstat(6,1));
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
ylabel('% \Delta Prod P')
title('min % \Delta Prod P')
axis([-0.5 0.3 -0.05 0.02])

subplot(2,2,4)
scatter(dNlP(:,np), dPlP(:,np), 10, ptemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('% \Delta Prod P')
title('min % \Delta Prod P')
axis([-0.3 0.5 -0.05 0.02])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_nucmax_vs_prod_minP_pdiff.png'])


%% LD max
xd = find(psets==pdstat(7,1));
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
ylabel('% \Delta Prod D')
title('max % \Delta Prod D')
axis([-0.05 0.3 -0.015 0.01])

subplot(2,2,4)
scatter(dNlD(:,xd), dPlD(:,xd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('% \Delta Prod D')
title('max % \Delta Prod D')
axis([-0.07 0.08 -0.015 0.01])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_nucmax_vs_prod_maxD_pdiff.png'])

%% LD min
nd = find(psets==pdstat(8,1));
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
ylabel('% \Delta Prod D')
title('min % \Delta Prod D')
axis([-0.05 0.3 -0.015 0.01])

subplot(2,2,4)
scatter(dNlD(:,nd), dPlD(:,nd), 10, btemp_fore, 'filled')
cmocean('thermal')
colorbar
caxis([0 30])
xlabel('\Delta Nu/Cmax')
ylabel('% \Delta Prod D')
title('min % \Delta Prod D')
axis([-0.07 0.08 -0.015 0.01])
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_nucmax_vs_prod_minD_pdiff.png'])

