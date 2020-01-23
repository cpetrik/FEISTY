% Time series of TEeff
% Need: AllL mdet  mmz_loss  mlz_loss mnpp
% Ensemble results

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% Zoop and det and npp 
load([gpath 'cobalt_det_temp_zoop_npp_means.mat']);

% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];


nfile = 'param_ensemble/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs';
ppath = [pp nfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mzloss_5yr_hist = mzloss_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_5yr_hist = lzloss_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_hist = det_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_hist = npp_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

hmmz_loss = mzloss_5yr_hist;
hmlz_loss = lzloss_5yr_hist;
hmdet = det_5yr_hist;
hmnpp = npp_5yr_hist;

mzloss_5yr_fore = mzloss_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzloss_5yr_fore = lzloss_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_fore = det_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_fore = npp_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

fmmz_loss = mzloss_5yr_fore;
fmlz_loss = lzloss_5yr_fore;
fmdet = det_5yr_fore;
fmnpp = npp_5yr_fore;

%% Hindcast
% 5yr means at each grid cell
load([fpath 'Means_Historic_',harv,'_prod_' cfile '.mat'],'lp_prod','ld_prod');
HL = lp_prod+ld_prod;
clear lp_prod ld_prod

%% RCP 8.5
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'lp_prod','ld_prod');
FL = lp_prod+ld_prod;
clear lp_prod ld_prod

%% Ensemble parameter sets
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([epath 'Historic_All_fish03_ensem6_mid_bestAIC_multFup_multPneg.mat'],...
    'htTEeff_LTL','htTEeff_HTL','htTEeff_ATL','htTE_HTL','htTE_ATL');
load([epath 'Forecast_All_fish03_ensem6_mid_bestAIC_multFup_multPneg.mat'],...
    'ftTEeff_LTL','ftTEeff_HTL','ftTEeff_ATL','ftTE_HTL','ftTE_ATL');

%% ts
y1 = 1860+(1/12):5:2005;
y2 = 2005+(1/12):5:2100;
y = [y1 y2];

L = [HL FL];
L(L<0) = 0;
npp = [hmnpp fmnpp];
det = [hmdet fmdet];
mz_loss = [hmmz_loss fmmz_loss];
lz_loss = [hmlz_loss fmlz_loss];

mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_loss = nanmean(mz_loss);
mlz_loss = nanmean(lz_loss);
mL = nanmean(L);

%% Effective TEs
% With BE*det instead of Bent
%TEeff_ATL = production_L/NPP
TEeff_ATL = L./npp;
TEeff_ATL(TEeff_ATL==-Inf) = NaN;
TEeff_ATL(TEeff_ATL==Inf) = NaN;
TEeff_ATL(TEeff_ATL<0) = NaN;
TEeff_ATL(TEeff_ATL>1) = NaN;
%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTL = (BE*det + mz_loss + lz_loss)./npp;
TEeff_LTL(TEeff_LTL==-Inf) = NaN;
TEeff_LTL(TEeff_LTL==Inf) = NaN;
TEeff_LTL(TEeff_LTL<0) = NaN;
TEeff_LTL(TEeff_LTL>1) = NaN;
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTL = L./(BE*det + mz_loss + lz_loss); 
TEeff_HTL(TEeff_HTL<0) = NaN;
TEeff_HTL(TEeff_HTL>1) = NaN;

%Calc TE from effTEs
TELTL = real(TEeff_LTL.^(1/(4/3)));
TEATL = real(TEeff_ATL.^(1/4));         
TEHTL = real(TEeff_HTL.^(1/3));

% Means
q(1,:) = nanmean(TEeff_LTL);
q(2,:) = nanmean(TEeff_HTL);
q(3,:) = nanmean(TEeff_ATL);
q(4,:) = nanmean(TELTL);
q(5,:) = nanmean(TEHTL);
q(6,:) = nanmean(TEATL);

%Ensemble
eTEeff_ATL = [htTEeff_ATL ftTEeff_ATL];
eTEeff_HTL = [htTEeff_HTL ftTEeff_HTL];
eTE_ATL = [htTE_ATL ftTE_ATL];
eTE_HTL = [htTE_HTL ftTE_HTL];

%% save
save([epath 'TEeffDet_5yr_means_Hist_Fore_All_fish03_ensem6_mid_bestAIC_multFup_multPneg.mat'],...
    'L','mdet','mmz_loss','mlz_loss','mnpp','y',...
    'TEeff_ATL','TEeff_LTL','TEeff_HTL','eTE_ATL','eTE_HTL',...
    'eTEeff_ATL','eTEeff_HTL');

%% subplot of TE effs
figure(1)
subplot(2,2,1)
plot(y,100*q(4,:),'k','Linewidth',2); hold on;
title('TEeff LTL')
ylim([9 10])
xlim([1900 2095])
xlabel('Year')

subplot(2,2,2)
plot(y,100*eTEeff_HTL); hold on;
plot(y,100*q(2,:),'k','Linewidth',2); hold on;
title('TEeff HTL')
%ylim([9 10])
xlim([1900 2095])
xlabel('Year')

subplot(2,2,3)
plot(y,100*eTEeff_ATL); hold on;
plot(y,100*q(3,:),'k','Linewidth',2); hold on;
title('TEeff ATL')
ylim([0.03 0.11])
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TEeff_sub_ensem.png'])

%% subplot of TEs
figure(2)
subplot(2,2,1)
plot(y,100*q(4,:),'k','Linewidth',2); hold on;
title('TE LTL')
ylim([9 10])
xlim([1900 2095])
xlabel('Year')

subplot(2,2,2)
plot(y,100*eTE_HTL); hold on;
plot(y,100*q(5,:),'k','Linewidth',2); hold on;
title('TE HTL')
xlim([1900 2095])
xlabel('Year')

subplot(2,2,3)
plot(y,100*eTE_ATL); hold on;
plot(y,100*q(6,:),'k','Linewidth',2); hold on;
title('TE ATL')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_sub_ensem.png'])

%% subplot of TEs
figure(3)
subplot(2,2,1)
plot(y,100*q(4,:),'color',[0.5 0 1],'Linewidth',2); hold on;
title('TE LTL')
ylim([9 10])
xlim([1900 2095])
xlabel('Year')

subplot(2,2,2)
plot(y,100*eTE_HTL,'color',[0 0.5 0.75]); hold on;
plot(y,100*q(5,:),'color',[0 0.5 0.75],'Linewidth',2); hold on;
title('TE HTL')
xlim([1900 2095])
xlabel('Year')

subplot(2,2,3)
plot(y,100*eTE_ATL,'k'); hold on;
plot(y,100*q(6,:),'k','Linewidth',2); hold on;
title('TE ATL')
xlim([1900 2095])
xlabel('Year')

subplot(2,2,4)
plot(y,100*q(4,:),'color',[0.5 0 1],'Linewidth',2); hold on;
plot(y,100*eTE_HTL,'color',[0 0.5 0.75]); hold on;
plot(y,100*q(5,:),'color',[0 0.5 0.75],'Linewidth',2); hold on;
plot(y,100*eTE_ATL,'k'); hold on;
plot(y,100*q(6,:),'k','Linewidth',2); hold on;
%legend('LTL','HTL','ATL')
%legend('location','northeast')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_sub2_ensem.png'])

%% CONE OF UNCERTAINTY
mHTL = mean([eTE_HTL; q(5,:)]);
mATL = mean([eTE_ATL; q(6,:)]);

%find sims with min and max values
hmax = find(eTE_HTL(:,1) == max(eTE_HTL(:,1)));
hmin = find(eTE_HTL(:,1) == min(eTE_HTL(:,1)));
amax = find(eTE_ATL(:,1) == max(eTE_ATL(:,1)));
amin = find(eTE_ATL(:,1) == min(eTE_ATL(:,1)));

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%replace with sims with min and max values
Ya=[eTE_ATL(amax,:) fliplr(eTE_ATL(amin,:))]; 
Yh=[eTE_HTL(hmax,:) fliplr(eTE_HTL(hmin,:))]; 
 
%% subplot of TEs with cone
figure(4)
subplot(2,2,1)
plot(y,100*q(4,:),'k','Linewidth',2); hold on;
title('TE LTL')
ylim([9 10])
xlim([1900 2095])
xlabel('Year')

subplot(2,2,2)
fill(X,100*Yh,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*mHTL,'k','Linewidth',2); hold on;
title('TE HTL')
xlim([1900 2095])
xlabel('Year')

subplot(2,2,3)
fill(X,100*Ya,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*mATL,'k','Linewidth',2); hold on;
title('TE ATL')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_sub_cone_ensem.png'])

%% subplot of TEs with cone, 4th all together
figure(5)
subplot(2,2,1)
plot(y,100*q(4,:),'color',[0.5 0 1],'Linewidth',2); hold on;
title('TE LTL')
ylim([9 10])
xlim([1900 2095])
xlabel('Year')

subplot(2,2,2)
fill(X,100*Yh,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*mHTL,'color',[0 0.5 0.75],'Linewidth',2); hold on;
title('TE HTL')
xlim([1900 2095])
xlabel('Year')

subplot(2,2,3)
fill(X,100*Ya,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*mATL,'k','Linewidth',2); hold on;
title('TE ATL')
xlim([1900 2095])
xlabel('Year')

subplot(2,2,4)
plot(y,100*q(4,:),'color',[0.5 0 1],'Linewidth',2); hold on;
fill(X,100*Yh,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*mHTL,'color',[0 0.5 0.75],'Linewidth',2); hold on;
fill(X,100*Ya,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
plot(y,100*mATL,'k','Linewidth',2); hold on;
%legend('LTL','HTL','ATL')
%legend('location','northeast')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_sub2_cone_ensem.png'])
