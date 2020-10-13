% Visualize time series output of FEISTY ensembles
% Historic (1860-2005) and Forecast time period (2006-2100) at all locations
% Saved as mat files
% Ensemble mid6, temp3
% Find ensemble end members with max and min change in production relative 
% to their 1951 value

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
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
load([fpath 'ESM2M_Hist_Fore/Time_Means_Historic_Forecast_',harv,'_' cfile '.mat'],...
    'HF_tamean','HP_tamean','HD_tamean','HB_tamean',...
    'FF_tamean','FP_tamean','FD_tamean','FB_tamean',...
    'HA_tamean','FA_tamean');

%% Ensemble parameter sets
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'haTsF','haTsP','haTsD','haTmF','haTmP','haTmD','haTB','haTlP','haTlD');
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'faTsF','faTsP','faTsD','faTmF','faTmP','faTmD','faTB','faTlP','faTlD');

%aT=tamean=nansum(prod.*area,1);

% lpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',...
%     'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext');

%% calculate: 
%NO benthic production = benthic biomass * detritus flux * benthic efficiency
%benthic production = detritus flux * benthic efficiency
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([bpath 'cobalt_det_ts.mat'],'det_mean_hist','det_mean_fore',...
    'mo_hist','mo_fore');
det_hist = det_mean_hist' * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_fore = det_mean_fore' * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

% HB = haTB .* repmat(det_hist,43,1) .* 0.075;
% FB = faTB .* repmat(det_fore,43,1) .* 0.075;
% HB_tam = HB_tamean .* det_hist .* 0.075;
% FB_tam = FB_tamean .* det_fore .* 0.075;

HB = repmat(det_hist,43,1) .* 0.075;
FB = repmat(det_fore,43,1) .* 0.075;
HB_tam = det_hist .* 0.075;
FB_tam = det_fore .* 0.075;

%% ts
%In original saved file
y1 = 1860+(1/12):(1/12):2005;
y2 = 2005+(1/12):(1/12):2100;
y = [y1 y2];

HF = haTsF + haTmF;
HP = haTsP + haTmP + haTlP;
HD = haTsD + haTmD + haTlD;
HA = HF + HP + HD;

FF = faTsF + faTmF;
FP = faTsP + faTmP + faTlP;
FD = faTsD + faTmD + faTlD;
FA = FF + FP + FD;

tForig = [HF_tamean FF_tamean];
tPorig = [HP_tamean FP_tamean];
tDorig = [HD_tamean FD_tamean];
tAorig = [HA_tamean FA_tamean];
tBorig = [HB_tam FB_tam];

tF = [HF FF];
tP = [HP FP];
tD = [HD FD];
tA = [HA FA];
tB = [HB FB];
%tB = [haTB faTB];

tF = [tF; tForig];
tP = [tP; tPorig];
tD = [tD; tDorig];
tA = [tA; tAorig];
tB = [tB; tBorig];

mtF = mean(tF);
mtP = mean(tP);
mtD = mean(tD);
mtA = mean(tA);
mtB = mean(tB);

%% Percent Difference
hyr = find(y>1950 & y<=2000);
fyr = find(y>2050 & y<=2100);

% percent difference from 1951?
pdF = (tF(:,hyr(1)) - tF(:,end)) ./ tF(:,hyr(1)); 
pdP = (tP(:,hyr(1)) - tP(:,end)) ./ tP(:,hyr(1));
pdD = (tD(:,hyr(1)) - tD(:,end)) ./ tD(:,hyr(1)); 
pdA = (tA(:,hyr(1)) - tA(:,end)) ./ tA(:,hyr(1));
pdB = (tB(:,hyr(1)) - tB(:,end)) ./ tB(:,hyr(1));


%% Individual Plots colored by pdiff All fish
red_params(44,:) = [0.7, 0.175, 0.2, 4, 70, 0.0855];
se = [0.02, 0.001, 0.007, 0.1, 1, 0.001];
rparam = red_params + repmat(se,44,1).*randn(size(red_params));
%define scttered k
rparam(:,7) = 0.063 + 0.001*randn(44,1);

cmR=cbrewer('seq','Reds',50,'PCHIP');

%% aM vs. bM
figure(1)
%F
subplot(2,2,1)
scatter(rparam(:,2),rparam(:,4),30,100*pdF,'filled','MarkerFaceAlpha',0.7)
%colormap('copper')
colormap(cmR)
colorbar
xlim([0.14 0.21])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Temp-dep of metabolism (b_M)')
title('% \Delta F')

%P
subplot(2,2,2)
scatter(rparam(:,2),rparam(:,4),30,100*pdP,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.14 0.21])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Temp-dep of metabolism (b_M)')
title('% \Delta P')

%All
subplot(2,2,3)
scatter(rparam(:,2),rparam(:,4),30,100*pdD,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.14 0.21])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Temp-dep of metabolism (b_M)')
title('% \Delta D')

%All
subplot(2,2,4)
scatter(rparam(:,2),rparam(:,4),30,100*pdA,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.14 0.21])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Temp-dep of metabolism (b_M)')
title('% \Delta All')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_aM_bM_prod_pdiff.png'])


%% aM vs. assim
figure(2)
%F
subplot(2,2,1)
scatter(rparam(:,1),rparam(:,4),30,100*pdF,'filled','MarkerFaceAlpha',0.7)
%colormap('copper')
colormap(cmR)
colorbar
xlim([0.55 0.8])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Assimilation effic (\alpha)')
title('% \Delta F')

%P
subplot(2,2,2)
scatter(rparam(:,1),rparam(:,4),30,100*pdP,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.55 0.8])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Assimilation effic (\alpha)')
title('% \Delta P')

%All
subplot(2,2,3)
scatter(rparam(:,1),rparam(:,4),30,100*pdD,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.55 0.8])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Assimilation effic (\alpha)')
title('% \Delta D')

%All
subplot(2,2,4)
scatter(rparam(:,1),rparam(:,4),30,100*pdA,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.55 0.8])
ylim([2.5 5.5])
ylabel('Weight-dep of metabolism (a_M)')
xlabel('Assimilation effic (\alpha)')
title('% \Delta All')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_aM_assim_prod_pdiff.png'])

%% k vs. kt
figure(3)
%F
subplot(2,2,1)
scatter(rparam(:,7),rparam(:,6),30,100*pdF,'filled','MarkerFaceAlpha',0.7)
%colormap('copper')
colormap(cmR)
colorbar
xlim([0.06 0.066])
ylim([0.07 0.1])
ylabel('Temp-dep metabolism (k_M)')
xlabel('Temp-dep feeding (k)')
title('% \Delta F')

%P
subplot(2,2,2)
scatter(rparam(:,7),rparam(:,6),30,100*pdP,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.06 0.066])
ylim([0.07 0.1])
ylabel('Temp-dep metabolism (k_M)')
xlabel('Temp-dep feeding (k)')
title('% \Delta P')

%All
subplot(2,2,3)
scatter(rparam(:,7),rparam(:,6),30,100*pdD,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.06 0.066])
ylim([0.07 0.1])
ylabel('Temp-dep metabolism (k_M)')
xlabel('Temp-dep feeding (k)')
title('% \Delta D')

%All
subplot(2,2,4)
scatter(rparam(:,7),rparam(:,6),30,100*pdA,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.06 0.066])
ylim([0.07 0.1])
ylabel('Temp-dep metabolism (k_M)')
xlabel('Temp-dep feeding (k)')
title('% \Delta All')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_k_kt_prod_pdiff.png'])

%% aE vs. bE
figure(4)
%F
subplot(2,2,1)
scatter(rparam(:,3),rparam(:,5),30,100*pdF,'filled','MarkerFaceAlpha',0.7)
%colormap('copper')
colormap(cmR)
colorbar
xlim([0.125 0.275])
ylim([40 110])
ylabel('Weight-dep of encounter (a_E)')
xlabel('Temp-dep of encounter (b_E)')
title('% \Delta F')

%P
subplot(2,2,2)
scatter(rparam(:,3),rparam(:,5),30,100*pdP,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.125 0.275])
ylim([40 110])
ylabel('Weight-dep of encounter (a_E)')
xlabel('Temp-dep of encounter (b_E)')
title('% \Delta P')

%All
subplot(2,2,3)
scatter(rparam(:,3),rparam(:,5),30,100*pdD,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.125 0.275])
ylim([40 110])
ylabel('Weight-dep of encounter (a_E)')
xlabel('Temp-dep of encounter (b_E)')
title('% \Delta D')

%All
subplot(2,2,4)
scatter(rparam(:,3),rparam(:,5),30,100*pdA,'filled','MarkerFaceAlpha',0.7)
colormap(cmR)
colorbar
xlim([0.125 0.275])
ylim([40 110])
ylabel('Weight-dep of encounter (a_E)')
xlabel('Temp-dep of encounter (b_E)')
title('% \Delta All')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ensem_mid6_temp3_aE_bE_prod_pdiff.png'])


