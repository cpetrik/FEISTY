%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';


%% POEM params
bees = 0.175:0.005:0.195;

kays = 0.0405:0.01:0.125;
kt = kays(5);
tkfn = num2str(100+int64(100*kt));

Fish = 0.1:0.1:0.6;
frate = Fish(3);
tfish = num2str(100+int64(10*frate));

charv = ['All_fish',tfish(2:end)];
harv = tfish(2:end);

load([dp 'Clim_k08_fished',harv,'_LME_SAUP_catch_comp_b.mat']);
rall8 = rall;
rF8 = rF;
rP8 = rP;
rD8 = rD;
rmse8 = rmse;
rmseF8 = rmseF;
rmseP8 = rmseP;
rmseD8 = rmseD;
clear rall rF rP rD rmse rmseF rmseP rmseD

load([dp 'Clim_k09_fished',harv,'_LME_SAUP_catch_comp_b.mat']);
rall9 = rall;
rF9 = rF;
rP9 = rP;
rD9 = rD;
rmse9 = rmse;
rmseF9 = rmseF;
rmseP9 = rmseP;
rmseD9 = rmseD;
clear rall rF rP rD rmse rmseF rmseP rmseD

%% Plots of r and RMSE
cfile2 = ['Clim_fished',harv,'_Dc_enc70_cmax-metab20_k0809_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];

figure(2)
subplot(2,2,1)
plot(bees,rall8,'r'); hold on;
plot(bees,rall9,'b');
legend('.08','.09')
xlim([bees(1) bees(end)])
ylabel('r')
xlabel('Temp exponent')
title('Corr all fish')

subplot(2,2,2)
plot(bees,rF8,'r'); hold on;
plot(bees,rF9,'b');
xlim([bees(1) bees(end)])
ylabel('r')
xlabel('Temp exponent')
title('Corr F')

subplot(2,2,3)
plot(bees,rP8,'r'); hold on;
plot(bees,rP9,'b');
xlim([bees(1) bees(end)])
ylabel('r')
xlabel('Temp exponent')
title('Corr P')

subplot(2,2,4)
plot(bees,rD8,'r'); hold on;
plot(bees,rD9,'b');
xlim([bees(1) bees(end)])
ylabel('r')
xlabel('Temp exponent')
title('Corr D')
stamp(cfile2)
print('-dpng',[pp 'Clim_k0809_fished',harv,'_LME_SAUP_catch_corr_b'])


%%

figure(3)
subplot(2,2,1)
plot(bees,rmse8,'r'); hold on;
plot(bees,rmse9,'b');
legend('.08','.09')
xlim([bees(1) bees(end)])
ylabel('RMSE')
xlabel('Temp exponent')
title('RMSE all fish')

subplot(2,2,2)
plot(bees,rmseF8,'r'); hold on;
plot(bees,rmseF9,'b');
xlim([bees(1) bees(end)])
ylabel('RMSE')
xlabel('Temp exponent')
title('RMSE F')

subplot(2,2,3)
plot(bees,rmseP8,'r'); hold on;
plot(bees,rmseP9,'b');
xlim([bees(1) bees(end)])
ylabel('RMSE')
xlabel('Temp exponent')
title('RMSE P')

subplot(2,2,4)
plot(bees,rmseD8,'r'); hold on;
plot(bees,rmseD9,'b');
xlim([bees(1) bees(end)])
ylabel('RMSE')
xlabel('Temp exponent')
title('RMSE D')
stamp(cfile2)
print('-dpng',[pp 'Clim_k0809_fished',harv,'_LME_SAUP_catch_rmse_b'])




