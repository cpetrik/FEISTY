% Check that new version with param structure gives same results

clear all
close all

%% Last year means
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath = ['/Volumes/FEISTY/NC/Clim_comp_tests/' cfile '/'];
load([fpath,'Orig_Means_Climatol_All_fish03_' cfile '.mat'])

%%
sp_mean1 = sp_mean;
sf_mean1 = sf_mean;
sd_mean1 = sd_mean;
mp_mean1 = mp_mean;
mf_mean1 = mf_mean;
md_mean1 = md_mean;
lp_mean1 = lp_mean;
ld_mean1 = ld_mean;
b_mean1  = b_mean;

%% Full climatol
%/Volumes/FEISTY/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/Climatology
fpath2 = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/Climatology/'];
load([fpath2,'Means_Climatol_All_fish03_' cfile '.mat'])

sp_mean2 = sp_mean;
sf_mean2 = sf_mean;
sd_mean2 = sd_mean;
mp_mean2 = mp_mean;
mf_mean2 = mf_mean;
md_mean2 = md_mean;
lp_mean2 = lp_mean;
ld_mean2 = ld_mean;
b_mean2  = b_mean;

%% Same as full climatol?
sfs = sum(sf_mean1==sf_mean2);
sps = sum(sp_mean1==sp_mean2);
sds = sum(sd_mean1==sd_mean2);
mfs = sum(mf_mean1==mf_mean2);
mps = sum(mp_mean1==mp_mean2);
mds = sum(md_mean1==md_mean2);
lps = sum(lp_mean1==lp_mean2);
lds = sum(ld_mean1==ld_mean2);
bs = sum(b_mean1==b_mean2);

%% Plot
x=0.1:0.5:50;
subplot(3,3,1)
plot(sf_mean2,sf_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 1.5 0 1.5])
title('SF')

subplot(3,3,2)
plot(sp_mean2,sp_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 0.5 0 0.5])
title('SP')

subplot(3,3,3)
plot(sd_mean2,sd_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 0.75 0 0.75])
title('SD')

subplot(3,3,4)
plot(mf_mean2,mf_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 25 0 25])
title('MF')

subplot(3,3,5)
plot(mp_mean2,mp_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 15 0 15])
title('MP')

subplot(3,3,6)
plot(md_mean2,md_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 10 0 10])
title('MD')

subplot(3,3,7)
plot(lp_mean2,lp_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 25 0 25])
title('LP')

subplot(3,3,8)
plot(ld_mean2,ld_mean1,'.k'); hold on;
plot(x,x,'b')
axis([0 50 0 50])
title('LD')

subplot(3,3,9)
plot(b_mean2,b_mean1,'.k'); hold on;
plot(x,x,'b')
%axis([0 80 0 50])
title('B')

%% Last year means 
fname3 = ['/Users/cpetrik/Dropbox/Princeton/POEM_other/FEISTY_clim_complete/',...
    'Output/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/Climatol_fish_F030_P010_D030'];

