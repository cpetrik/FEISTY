% Time series of TEeff
% Need: AllL mdet  mmz_prod  mlz_prod mnpp
% Ensemble results mid6, same k's

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

% Zoop and det and npp 
load([gpath 'cobalt_det_temp_zoop_npp_means.mat'],...
    'ptemp_1yr_hist','ptemp_1yr_fore',...
    'btemp_1yr_hist','btemp_1yr_fore',...
    'det_1yr_hist','det_1yr_fore',...
    'npp_1yr_hist','npp_1yr_fore',...
    'lzprod_1yr_hist','lzprod_1yr_fore',...
    'mzprod_1yr_hist','mzprod_1yr_fore');

% FEISTY
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
efile = 'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';
sfile = 'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',efile,'/full_runs/'];
%Orig
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
fpath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];

% Ensemble parameter sets w/baseline
epath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',efile,'/'];
epath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',efile,'/'];
spath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',sfile,'/'];
spath2 = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/',sfile,'/'];

BE = 0.075;

%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mzprod_1yr_hist = mzprod_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_1yr_hist = lzprod_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_1yr_hist = det_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_1yr_hist = npp_1yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

hmmz_prod = mzprod_1yr_hist;
hmlz_prod = lzprod_1yr_hist;
hmdet = det_1yr_hist;
hmnpp = npp_1yr_hist;

mzprod_1yr_fore = mzprod_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_1yr_fore = lzprod_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_1yr_fore = det_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_1yr_fore = npp_1yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

fmmz_prod = mzprod_1yr_fore;
fmlz_prod = lzprod_1yr_fore;
fmdet = det_1yr_fore;
fmnpp = npp_1yr_fore;

%% Hindcast
% 1yr means at each grid cell
load([fpath 'Historic_ESM2M/Means_Historic_',harv,'_prod_' cfile '.mat'],'lp_prod1','ld_prod1');
HL = lp_prod1+ld_prod1;
clear lp_prod1 ld_prod1

%% RCP 8.5
load([fpath 'Forecast_RCP85_ESM2M/Means_fore_',harv,'_' cfile '.mat'],'lp_prod1','ld_prod1');
FL = lp_prod1+ld_prod1;
clear lp_prod1 ld_prod1

%% Ensemble parameter sets
% Vary temp-dep
load([epath 'Historic_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'ht1TEeff_HTL','ht1TEeff_ATL','ht1TE_HTL','ht1TE_ATL');
load([epath 'Forecast_All_fish03_ensem6_mid_temp3_bestAIC_multFup_multPneg.mat'],...
    'ft1TEeff_HTL','ft1TEeff_ATL','ft1TE_HTL','ft1TE_ATL');
eTEeff_ATL = [ht1TEeff_ATL ft1TEeff_ATL];
eTEeff_HTL = [ht1TEeff_HTL ft1TEeff_HTL];
eTE_ATL = [ht1TE_ATL ft1TE_ATL];
eTE_HTL = [ht1TE_HTL ft1TE_HTL];

% Equal temp-dep
load([spath 'Historic_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'ht1TEeff_LTL','ht1TEeff_HTL','ht1TEeff_ATL','ht1TE_HTL','ht1TE_ATL');
load([spath 'Forecast_All_fish03_ensem6_mid_samek_bestAIC_multFup_multPneg.mat'],...
    'ft1TEeff_LTL','ft1TEeff_HTL','ft1TEeff_ATL','ft1TE_HTL','ft1TE_ATL');
sTEeff_ATL = [ht1TEeff_ATL ft1TEeff_ATL];
sTEeff_HTL = [ht1TEeff_HTL ft1TEeff_HTL];
sTE_ATL = [ht1TE_ATL ft1TE_ATL];
sTE_HTL = [ht1TE_HTL ft1TE_HTL];

%% ts
y1 = 1860+(1/12):1:2005;
y2 = 2005+(1/12):1:2100;
y = [y1 y2];

L = [HL FL];
L(L<0) = 0;
npp = [hmnpp fmnpp];
det = [hmdet fmdet];
mz_prod = [hmmz_prod fmmz_prod];
lz_prod = [hmlz_prod fmlz_prod];

mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_prod = nanmean(mz_prod);
mlz_prod = nanmean(lz_prod);
mL = nanmean(L);

%% Effective TEs
% With BE*det instead of Bent
% With Zprod instead of Zloss
%TEeff_ATL = production_L/NPP
TEeff_ATL = L./npp;
TEeff_ATL(TEeff_ATL==-Inf) = NaN;
TEeff_ATL(TEeff_ATL==Inf) = NaN;
TEeff_ATL(TEeff_ATL<0) = NaN;
TEeff_ATL(TEeff_ATL>1) = NaN;
%TEeff_LTL = (production_benthic_invert+mesozoo_prod)/NPP
TEeff_LTL = (BE*det + mz_prod + lz_prod)./npp;
TEeff_LTL(TEeff_LTL==-Inf) = NaN;
TEeff_LTL(TEeff_LTL==Inf) = NaN;
TEeff_LTL(TEeff_LTL<0) = NaN;
TEeff_LTL(TEeff_LTL>1) = NaN;
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod)
TEeff_HTL = L./(BE*det + mz_prod + lz_prod); 
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


%% save
save([epath 'TEeff_Det_Zprod_1yr_means_Hist_Fore_All_fish03_ensem6_mid_temp3_samek_bestAIC_multFup_multPneg.mat'],...
    'L','mdet','mmz_prod','mlz_prod','mnpp','y',...
    'TEeff_ATL','TEeff_LTL','TEeff_HTL',...
    'eTE_ATL','eTE_HTL','eTEeff_ATL','eTEeff_HTL',...
    'sTE_ATL','sTE_HTL','sTEeff_ATL','sTEeff_HTL');
save([epath2 'TEeff_Det_Zprod_1yr_means_Hist_Fore_All_fish03_ensem6_mid_temp3_samek_bestAIC_multFup_multPneg.mat'],...
    'L','mdet','mmz_prod','mlz_prod','mnpp','y',...
    'TEeff_ATL','TEeff_LTL','TEeff_HTL',...
    'eTE_ATL','eTE_HTL','eTEeff_ATL','eTEeff_HTL',...
    'sTE_ATL','sTE_HTL','sTEeff_ATL','sTEeff_HTL');

% load([epath 'TEeff_Det_Zprod_1yr_means_Hist_Fore_All_fish03_ensem6_mid_temp3_samek_bestAIC_multFup_multPneg.mat']);

%% CONE OF UNCERTAINTY
% Add to Ensemble
%q(4,:) = nanmean(TELTL);
teLTL = q(4,:);
teHTL = [eTE_HTL; q(5,:)];
teATL = [eTE_ATL; q(6,:)];

% difference from 1951
id=find(y>1950);
yid=id(1);
dLTL = teLTL - teLTL(yid);
deHTL = teHTL - teHTL(:,yid);
deATL = teATL - teATL(:,yid);
dsHTL = sTE_HTL - sTE_HTL(:,yid);
dsATL = sTE_ATL - sTE_ATL(:,yid);

meHTL = mean(deHTL);
meATL = mean(deATL);
seHTL = std(deHTL);
seATL = std(deATL);

msHTL = mean(dsHTL);
msATL = mean(dsATL);
ssHTL = std(dsHTL);
ssATL = std(dsATL);

%create continuous x value array for plotting
X=[y fliplr(y)]; 
%create y values for out and then back
%+/- 1 stdev
eSa=[meATL+seATL fliplr(meATL-seATL)]; 
eSh=[meHTL+seHTL fliplr(meHTL-seHTL)];

sSa=[msATL+ssATL fliplr(msATL-ssATL)]; 
sSh=[msHTL+ssHTL fliplr(msHTL-ssHTL)];
 

%% subplot of TEs with cone stdev all together color
figure(1)
subplot(2,2,1)
plot(y,100*dLTL,'color',[0.5 0 1],'Linewidth',2); hold on;
plot(y,100*meHTL,'color',[0 0.5 0.75],'Linewidth',2); hold on;
plot(y,100*meATL,'k','Linewidth',2); hold on;
fill(X,100*eSh,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*eSa,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
ylim([-3.5 1.5])
xlim([1951 2095])
xlabel('Year')
ylabel('Percent change in TE relative to 1951')
text(1951,1.7,'A')

subplot(2,2,2)
plot(y,100*dLTL,'color',[0.5 0 1],'Linewidth',2); hold on;
plot(y,100*msHTL,'color',[0 0.5 0.75],'Linewidth',2); hold on;
plot(y,100*msATL,'k','Linewidth',2); hold on;
fill(X,100*sSh,'c','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
fill(X,100*sSa,'k','FaceAlpha',0.25,'EdgeAlpha',0.25); hold on; %plot filled area
legend('LTL','HTL','ATL')
legend('location','southwest')
ylim([-3.5 1.5])
xlim([1951 2095])
xlabel('Year')
text(1951,1.7,'B')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_Det_Zprod_d1951_cone_1std_ensem_both_color.png'])


%%
figure(2)
for i=1:15
    subplot(5,3,i)
    plot(y,sTE_HTL(i,:))
    ylim([0.13 0.23])
end

%% single plot of TEs with baseline for ppt
figure(3)
plot(y,100*dLTL,'color',[0.5 0 1],'Linewidth',3); hold on;
plot(y,100*meHTL,'color',[0 0.5 0.75],'Linewidth',3); hold on;
plot(y,100*meATL,'k','Linewidth',3); hold on;
ylim([-3 0.5])
xlim([1951 2095])
xlabel('Year')
ylabel('Percent change in TE relative to 1951')
ax=gca;
ax.FontSize = 16;
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_Det_Zprod_d1951_solo_color.png'])

