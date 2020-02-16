% Time series of TEeff
% Use Det & Zprod

clear all
close all

Pdir = '/Volumes/FEISTY/POEM_JLD/esm26_hist/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
gpath='/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); 
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
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

cmYOR=cbrewer('seq','YlOrRd',50,'PCHIP');
cmRP=cbrewer('seq','RdPu',50,'PCHIP');
cmPR=cbrewer('seq','PuRd',50,'PCHIP');


%% Zoop and det and npp 

%ESM2M in mmol N m-2 or mmol N m-2 d-1
% molN/m2/s --> g/m2/d
% 106/16 mol C in 1 mol N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W
mzprod_5yr_hist = mzprod_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_5yr_hist = lzprod_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_hist = det_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_hist = npp_5yr_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

hmmz_prod = mzprod_5yr_hist;
hmlz_prod = lzprod_5yr_hist;
hmdet = det_5yr_hist;
hmnpp = npp_5yr_hist;

mzprod_5yr_fore = mzprod_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_5yr_fore = lzprod_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
det_5yr_fore = det_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
npp_5yr_fore = npp_5yr_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

fmmz_prod = mzprod_5yr_fore;
fmlz_prod = lzprod_5yr_fore;
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
mz_prod = [hmmz_prod fmmz_prod];
lz_prod = [hmlz_prod fmlz_prod];


%% Effective TEs
% With BE*det instead of Bent
%TEeff_L = production_L/NPP
TEeff_L = L./npp;
TEeff_L(TEeff_L==-Inf) = NaN;
TEeff_L(TEeff_L==Inf) = NaN;
TEeff_L(TEeff_L<0) = NaN;
TEeff_L(TEeff_L>1) = NaN;
%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTL = (BE*det + mz_prod + lz_prod)./npp;
TEeff_LTL(TEeff_LTL==-Inf) = NaN;
TEeff_LTL(TEeff_LTL==Inf) = NaN;
TEeff_LTL(TEeff_LTL<0) = NaN;
TEeff_LTL(TEeff_LTL>1) = NaN;
%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTL = L./(BE*det + mz_prod + lz_prod); 
TEeff_HTL(TEeff_HTL<0) = NaN;
TEeff_HTL(TEeff_HTL>1) = NaN;

TELTL = real(TEeff_LTL.^(1/(4/3)));
TEL = real(TEeff_L.^(1/4));         
TEHTL = real(TEeff_HTL.^(1/3));   

%%

mnpp = nanmean(npp);
mdet = nanmean(det);
mmz_prod = nanmean(mz_prod);
mlz_prod = nanmean(lz_prod);
mL = nanmean(L);

cob(1,:) = nanmean(npp);
cob(2,:) = nanmean(det);
cob(3,:) = nanmean(mz_prod);
cob(4,:) = nanmean(lz_prod);
cob(5,:) = nanmean(L);

q(1,:) = nanmean(TEeff_LTL);
q(2,:) = nanmean(TEeff_HTL);
q(3,:) = nanmean(TEeff_L);
q(4,:) = nanmean(TELTL);
q(5,:) = nanmean(TEHTL);
q(6,:) = nanmean(TEL);


Q = array2table(q,'RowNames',{'TEeff_LTL','TEeff_HTL',...
    'TEeff_L','TE_LTL','TE_HTL','TE_ATL'});

%% save
writetable(Q,[fpath 'TEeffDetZprod_5yr_means_Hist_Fore_All_fish03_' cfile '.csv'],'Delimiter',',');

save([fpath 'TEeffDetZprod_5yr_means_Hist_Fore_All_fish03_' cfile '.mat'],...
    'L','mdet','mmz_prod','mlz_prod','mnpp','y',...
    'TEeff_L','TEeff_LTL','TEeff_HTL','q','Q','cob');

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
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_Zprod.png'])

%% subplot of TEs
figure(3)
subplot(2,2,1)
%plot(y,100*q(1,:),'k','Linewidth',2); hold on; %TEeff
plot(y,100*q(4,:),'k','Linewidth',2); hold on; %TE
title('TE LTL')
ylim([14 15.5])
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
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_Zprod_sub3.png'])

%% subplot of TEs
figure(4)
subplot(2,2,1)
plot(y,100*q(4,:),'color',[0.5 0 1],'Linewidth',2); hold on;
title('TE LTL')
ylim([14 15.5])
xlim([1900 2095])
xlabel('Year')
subplot(2,2,2)
plot(y,100*q(5,:),'color',[0 0.5 0.75],'Linewidth',2); hold on;
title('TE HTL')
xlim([1900 2095])
xlabel('Year')
subplot(2,2,3)
plot(y,100*q(6,:),'k','Linewidth',2); hold on;
title('TE ATL')
xlim([1900 2095])
xlabel('Year')
subplot(2,2,4)
plot(y,100*q(4,:),'color',[0.5 0 1],'Linewidth',2); hold on;
plot(y,100*q(5,:),'color',[0 0.5 0.75],'Linewidth',2); hold on;
plot(y,100*q(6,:),'k','Linewidth',2); hold on;
%legend('LTL','HTL','ATL')
%legend('location','northeast')
xlim([1900 2095])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TE_Zprod_sub4.png'])

%% subplot of TEeffs
figure(5)
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
print('-dpng',[ppath 'Hist_Fore_',harv,'_TEeff_Zprod_sub.png'])

%% subplot of TEeffs w/color
figure(2)
subplot(2,2,1)
plot(y,q(1,:),'color',[0.5 0 1],'Linewidth',2); hold on;
title('TEeff LTL')
xlim([1900 2095])
ylim([7.5e-2 8.5e-2])
xlabel('Year')
subplot(2,2,2)
plot(y,q(2,:),'color',[0 0.5 0.75],'Linewidth',2); hold on;
title('TEeff HTL')
xlim([1900 2095])
ylim([4.5e-3 7e-3])
xlabel('Year')
subplot(2,2,3)
plot(y,q(3,:),'k','Linewidth',2); hold on;
title('TEeff ATL')
xlim([1900 2095])
ylim([4.5e-4 7e-4])
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_TEeff_Zprod_sub_color.png'])

%% Prey to NPP
figure(6)
plot(y,(cob(2,:)./cob(1,:)),'k','Linewidth',2); hold on;
plot(y,(cob(3,:)./cob(1,:)),'r','Linewidth',2); hold on;
plot(y,(cob(4,:)./cob(1,:)),'b','Linewidth',2); hold on;
xlim([1900 2095])
title('Production / NPP')
legend('Det','MZ','LZ')
legend('location','west')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_Zprod_NPP_ratio.png'])

%% Zoop to Det
mz = mz_prod + lz_prod;
ZD = nanmean(mz ./ det);

figure(7)
plot(y,ZD,'k','Linewidth',2); hold on;
xlim([1900 2095])
title('Z / Det')
xlabel('Year')
print('-dpng',[ppath 'Hist_Fore_',harv,'_ZpDet_ratio.png'])


%%
figure(8)
subplot(2,2,1)
hist(log10(npp(:)))
title('NPP')
subplot(2,2,2)
hist(log10(det(:)))
title('Det')
subplot(2,2,3)
hist(log10(mz_prod(:) + lz_prod(:)))
title('Zoop')
subplot(2,2,4)
hist(log10(L(:)))
title('Large')

figure(9)
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





