% Time series of TEeff
% Need: AllL mdet  mmz_loss  mlz_loss mnpp

clear all
close all

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
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
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

cmYOR=cbrewer('seq','YlOrRd',50);
cmRP=cbrewer('seq','RdPu',50);
cmPR=cbrewer('seq','PuRd',50);


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
load([fpath 'Means_Historic_',harv,'_prod_' cfile '.mat'],'lp_prod','ld_prod');
HL = lp_prod+ld_prod;
clear lp_prod ld_prod

%% RCP 8.5
load([fpath 'Means_fore_',harv,'_' cfile '.mat'],'lp_prod','ld_prod');
FL = lp_prod+ld_prod;
clear lp_prod ld_prod

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


%% Effective TEs
% With BE*det instead of Bent
%TEeff_L = production_L/NPP
TEeff_L = L./npp;
TEeff_L(TEeff_L==-Inf) = NaN;
TEeff_L(TEeff_L==Inf) = NaN;
TEeff_L(TEeff_L<0) = NaN;
TEeff_L(TEeff_L>1) = NaN;
%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTLd = (BE*det + mz_loss + lz_loss)./npp;
TEeff_LTLd(TEeff_LTLd==-Inf) = NaN;
TEeff_LTLd(TEeff_LTLd==Inf) = NaN;
TEeff_LTLd(TEeff_LTLd<0) = NaN;
TEeff_LTLd(TEeff_LTLd>1) = NaN;
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTLd = L./(BE*det + mz_loss + lz_loss); 
TEeff_HTLd(TEeff_HTLd<0) = NaN;
TEeff_HTLd(TEeff_HTLd>1) = NaN;

TELTLd = real(TEeff_LTLd.^(1/(4/3)));
TEL = real(TEeff_L.^(1/4));         
TEHTLd = real(TEeff_HTLd.^(1/3));   

%%

mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_loss = nanmean(mz_loss);
mlz_loss = nanmean(lz_loss);
mL = nanmean(L);

cob(1,:) = nanmean(npp);
cob(2,:) = nanmean(det);
cob(3,:) = nanmean(mz_loss);
cob(4,:) = nanmean(lz_loss);
cob(5,:) = nanmean(L);

q(1,:) = nanmean(TEeff_LTLd);
q(2,:) = nanmean(TEeff_HTLd);
q(3,:) = nanmean(TEeff_L);
q(4,:) = nanmean(TELTLd);
q(5,:) = nanmean(TEHTLd);
q(6,:) = nanmean(TEL);


Q = array2table(q,'RowNames',{'TEeff_LTL','TEeff_HTL',...
    'TEeff_L','TE_LTL','TE_HTL','TE_ATL'});

%% save
writetable(Q,[fpath 'TEeffDet_5yr_means_Hist_Fore_All_fish03_' cfile '.csv'],'Delimiter',',');

save([fpath 'TEeffDet_5yr_means_Hist_Fore_All_fish03_' cfile '.mat'],...
    'L','mdet','mmz_loss','mlz_loss','mnpp','y',...
    'TEeff_L','TEeff_LTLd','TEeff_HTLd','q','Q','cob');

%%
figure(1)
plot(y,q(4,:),'k','Linewidth',2); hold on;
plot(y,q(5,:),'r','Linewidth',2); hold on;
plot(y,q(6,:),'b','Linewidth',2); hold on;
legend('LTL','HTL','ATL')
legend('location','northeast')
xlim([1900 2095])
xlabel('Year')
%ylabel('TE')
title(['TE'])
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE.png'])

%%
figure(2)
plot(y,log10(q(4,:)),'k','Linewidth',2); hold on;
plot(y,log10(q(5,:)),'r','Linewidth',2); hold on;
plot(y,log10(q(6,:)),'b','Linewidth',2); hold on;
legend('LTL','HTL','ATL')
legend('location','northeast')
xlim([1900 2095])
xlabel('Year')
%ylabel('TE')
title('log_1_0 TE')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_log10.png'])

%% subplot of TEs
figure(3)
subplot(2,2,1)
plot(y,100*q(4,:),'k','Linewidth',2); hold on;
title('TE LTL')
ylim([9 10])
xlim([1900 2095])
xlabel('Year')
subplot(2,2,2)
plot(y,100*q(5,:),'k','Linewidth',2); hold on;
title('TE HTL')
xlim([1900 2095])
xlabel('Year')
subplot(2,2,3)
plot(y,100*q(6,:),'k','Linewidth',2); hold on;
title('TE ATL')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_sub.png'])

%% subplot of TEeffs
figure(4)
subplot(2,2,1)
plot(y,q(1,:),'k','Linewidth',2); hold on;
title('TEeff LTL')
xlim([1900 2095])
xlabel('Year')
subplot(2,2,2)
plot(y,q(2,:),'k','Linewidth',2); hold on;
title('TEeff HTL')
xlim([1900 2095])
xlabel('Year')
subplot(2,2,3)
plot(y,q(3,:),'k','Linewidth',2); hold on;
title('TEeff ATL')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TEeff_sub.png'])

%% Prey to NPP
figure(5)
plot(y,(cob(2,:)./cob(1,:)),'k','Linewidth',2); hold on;
plot(y,(cob(3,:)./cob(1,:)),'r','Linewidth',2); hold on;
plot(y,(cob(4,:)./cob(1,:)),'b','Linewidth',2); hold on;
xlim([1900 2095])
title('Production / NPP')
legend('Det','MZ','LZ')
legend('location','southwest')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_prod_NPP_ratio.png'])

%% Zoop to Det
mz = mz_loss + lz_loss;
ZD = nanmean(mz ./ det);

figure(6)
plot(y,ZD,'k','Linewidth',2); hold on;
xlim([1900 2095])
title('Z / Det')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ZlDet_ratio.png'])


%%
figure(10)
subplot(2,2,1)
hist(log10(npp(:)))
title('NPP')
subplot(2,2,2)
hist(log10(det(:)))
title('Det')
subplot(2,2,3)
hist(log10(mz_loss(:) + lz_loss(:)))
title('Zoop')
subplot(2,2,4)
hist(log10(L(:)))
title('Large')

figure(11)
subplot(2,2,1)
plot(y,log10(cob(1,:)))
xlim([1900 2100])
title('NPP')
subplot(2,2,2)
plot(y,log10(cob(2,:)))
xlim([1900 2100])
title('Det')
subplot(2,2,3)
plot(y,log10(cob(3,:)),'b'); hold on;
plot(y,log10(cob(4,:)),'k'); hold on;
xlim([1900 2100])
title('Zoop')
subplot(2,2,4)
plot(y,log10(cob(5,:)))
xlim([1900 2100])
title('Large')
%%
figure(12)
subplot(2,2,1)
plot(y,(cob(2,:)./cob(1,:)))
xlim([1900 2100])
title('Det / NPP')
subplot(2,2,2)
plot(y,(cob(3,:)./cob(1,:)))
xlim([1900 2100])
title('M Zoop / NPP')
subplot(2,2,3)
plot(y,(cob(4,:)./cob(1,:)))
xlim([1900 2100])
title('L Zoop / NPP')
subplot(2,2,4)
plot(y,(cob(5,:)./cob(1,:)))
xlim([1900 2100])
title('Large / NPP')


%%
figure(20)
pcolor(mdet./mnpp)
shading flat
colorbar
caxis([-1 1])
cmocean('balance')

figure(21)
pcolor(mmz_loss./mnpp)
shading flat
colorbar
caxis([-1 1])
cmocean('balance')

figure(22)
pcolor(mlz_loss./mnpp)
shading flat
colorbar
caxis([-1 1])
cmocean('balance')

figure(23)
pcolor(L./mnpp)
shading flat
colorbar
caxis([-1 1])
cmocean('balance')



