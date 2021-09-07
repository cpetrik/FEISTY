% Visualize biomass of F and P as fn of zoop
% at 2 different temps (or bins)
% ESM2.6 Climatology of baseline params (kt=0.0855, k=0.063)
% and same k's (kt=0.063, k=0.063)

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Climatology grid
Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%% FEISTY baseline Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];
pp = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/' cfile '/'];

harv = 'All_fish03';

load([fpath 'Means_Climatol_',harv,'_' cfile '.mat'],...
        'sf_mean','sp_mean',...
        'mf_mean','mp_mean',...
        'lp_mean');

vF = sf_mean+mf_mean;
vP = sp_mean+mp_mean+lp_mean;

%% FEISTY equal k Output
efile = 'Dc_enc70-b200-k063_m4-b175-k063_c20-b250-k063_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100_kappaA50';
epath=['/Volumes/FEISTY/NC/Matlab_new_size/' efile '/'];

load([epath 'Means_Climatol_',harv,'_' efile '.mat'],...
        'sf_mean','sp_mean',...
        'mf_mean','mp_mean',...
        'lp_mean');

eF = sf_mean+mf_mean;
eP = sp_mean+mp_mean+lp_mean;

%% Zoop loss (don't have production)
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzloss_mean_clim',...
    'lzloss_mean_clim','ptemp_mean_clim');

% mgC/m2/d --> g/m2/d
mzprod_hist = mzloss_mean_clim * 1e-3 * 9.0;
lzprod_hist = lzloss_mean_clim * 1e-3 * 9.0;

z_hist = mzprod_hist + lzprod_hist;

%% vectorize
vT = ptemp_mean_clim(ID);
vZ = z_hist(ID);

%% mean temps around 10 and 25 C
lid = find(vT<11 & vT>=9);
hid = find(vT<26 & vT>=24);

%% fit lines varying
vZ10 = vZ(lid);
vZ20 = vZ(hid);
vF10 = vF(lid);
vF20 = vF(hid);
vP10 = vP(lid);
vP20 = vP(hid);

mF10 = fitlm(vZ10,vF10,'linear');
mP10 = fitlm(vZ10,vP10,'linear');
mF20 = fitlm(vZ20,vF20,'linear');
mP20 = fitlm(vZ20,vP20,'linear');

zoo = [0:0.25:3]';
pF10 = predict(mF10,zoo);
pP10 = predict(mP10,zoo);
pF20 = predict(mF20,zoo);
pP20 = predict(mP20,zoo);

opred(:,1) = pF10;
opred(:,2) = pP10;
opred(:,3) = pF20;
opred(:,4) = pP20;

%% fit lines equal
eF10 = eF(lid);
eF20 = eF(hid);
eP10 = eP(lid);
eP20 = eP(hid);

mF10 = fitlm(vZ10,eF10,'linear');
mP10 = fitlm(vZ10,eP10,'linear');
mF20 = fitlm(vZ20,eF20,'linear');
mP20 = fitlm(vZ20,eP20,'linear');

pF10 = predict(mF10,zoo);
pP10 = predict(mP10,zoo);
pF20 = predict(mF20,zoo);
pP20 = predict(mP20,zoo);

epred(:,1) = pF10;
epred(:,2) = pP10;
epred(:,3) = pF20;
epred(:,4) = pP20;

save([fpath 'FPbiomass_zProd_temp2_clim_fits_' cfile '.mat'],'opred','epred');

%%
figure(1)
subplot(2,2,1)
plot(vZ10,eP10,'o','color',[0/255 240/255 240/255]); hold on;
plot(vZ10,eF10,'x','color',[255/255 192/255 203/255]); hold on;
plot(zoo,epred(:,2),'b','LineWidth',2); hold on;
plot(zoo,epred(:,1),'r','LineWidth',2); hold on;
axis([0 3 0 20])
ylabel('Equal temp-dep')
title('10^oC')

subplot(2,2,2)
plot(vZ20,eP20,'o','color',[0/255 240/255 240/255]); hold on;
plot(vZ20,eF20,'x','color',[255/255 192/255 203/255]); hold on;
plot(zoo,epred(:,4),'b','LineWidth',2); hold on;
plot(zoo,epred(:,3),'r','LineWidth',2); hold on;
axis([0 3 0 20])
title('25^oC')

subplot(2,2,3)
plot(vZ10,vP10,'o','color',[0/255 240/255 240/255]); hold on;
plot(vZ10,vF10,'x','color',[255/255 192/255 203/255]); hold on;
plot(zoo,opred(:,2),'b','LineWidth',2); hold on;
plot(zoo,opred(:,1),'r','LineWidth',2); hold on;
ylabel('Unequal temp-dep')
axis([0 3 0 20])

%
subplot(2,2,4)
plot(vZ20,vP20,'o','color',[0/255 240/255 240/255]); hold on;
plot(vZ20,vF20,'x','color',[255/255 192/255 203/255]); hold on;
plot(zoo,opred(:,4),'b','LineWidth',2); hold on;
plot(zoo,opred(:,3),'r','LineWidth',2); hold on;
legend('P','F')
legend('location','northwest')
axis([0 3 0 20])
text(-0.5,-4,'Zooplankton production (g m^-^2 d^-^1)','HorizontalAlignment','center','FontWeight','bold')
text(-4.75,23,'Fish biomass (g m^-^2)','HorizontalAlignment','center','Rotation',90,'FontWeight','bold')
print('-dpng',[pp 'Clim_',harv,'_FPbiomass_zProd_temp2_fits.png'])




