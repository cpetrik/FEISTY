% P:D ratio by MEOWs
% Climatology
% 150 years
% Saved as mat files
% Compare to Daniel's model results & Reg Watson

clear all
close all

vpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

%% POEM
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';
ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%% DvD and Reg
load([vpath 'Ecoregion_comparison_v2_k85_BE075.mat'])

%% Comparison stats
nn1 = find(~isnan(f_large));
nn2 = find(~isnan(f_large_fwmodel));
nn3 = find(~isnan(f_catch));
nn0 = intersect(nn1,nn2);
nn = intersect(nn0,nn3);

f_large = f_large(nn);
f_large_fwmodel = f_large_fwmodel(nn);
f_catch = f_catch(nn);
Col_frac = Col_frac(nn);
Ecoregion = Ecoregion(nn);

diffD = Col_frac - f_large_fwmodel;
diffR = Col_frac - f_large;

id4 = find(f_catch<0.4);
id2 = find(f_catch<0.25);

%r
[rDall,pD]=corr(f_large_fwmodel,Col_frac);
[rDall2,pD2]=corr(f_large_fwmodel(id2),Col_frac(id2));
[rDall4,pD4]=corr(f_large_fwmodel(id4),Col_frac(id4));
[rRall,pR]=corr(f_large,Col_frac);
[rRall2,pR2]=corr(f_large(id2),Col_frac(id2));
[rRall4,pR4]=corr(f_large(id4),Col_frac(id4));

%root mean square error
% o=FracLP(did);
% p=plme_rPDcatch(did);
% n = length(o);
% num=nansum((p-o).^2);
% rmse = sqrt(num/n);


% Table
fish_stat(1,1) = rDall;
fish_stat(2,1) = rDall4;
fish_stat(3,1) = rDall2;
fish_stat(1,2) = rRall;
fish_stat(2,2) = rRall4;
fish_stat(3,2) = rRall2;

Fstat = array2table(fish_stat,'RowNames',{'All','0.40','0.25'},...
    'VariableNames',{'DvD','Reg'});
writetable(Fstat,[dpath 'Ecoreg_DvD_Reg_stats_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)
save([dpath 'Ecoreg_DvD_Reg_stats_' cfile '.mat'],'fish_stat')

%% Plot info
revamp=colormap(cmocean('amp'));
revamp=flipud(revamp);
close all

x=0:0.1:1;

%% Subplot with maps and corr
figure(1)
%DvD corr
subplot(3,3,1)
plot(f_large_fwmodel,Col_frac,'.k','MarkerSize',15); hold on;
plot(x,x,'--r');hold on;
text(0.05,1,['r = ' sprintf('%2.2f',rDall)])
axis([0 1.1 0 1.1])
xlabel('DvD')
ylabel('POEM')
title('All Ecoregions')

subplot(3,3,2)
plot(f_large_fwmodel(id4),Col_frac(id4),'.k','MarkerSize',15); hold on;
plot(x,x,'--r');hold on;
text(0.05,1,['r = ' sprintf('%2.2f',rDall4)])
axis([0 1.1 0 1.1])
xlabel('DvD')
ylabel('POEM')
title('UID < 0.40')

subplot(3,3,3)
plot(f_large_fwmodel(id2),Col_frac(id2),'.k','MarkerSize',15); hold on;
plot(x,x,'--r');hold on;
text(0.05,1,['r = ' sprintf('%2.2f',rDall2)])
axis([0 1.1 0 1.1])
xlabel('DvD')
ylabel('POEM')
title('UID < 0.25')

%Reg Corr
subplot(3,3,4)
plot(f_large,Col_frac,'.k','MarkerSize',15); hold on;
plot(x,x,'--r');hold on;
text(0.05,1,['r = ' sprintf('%2.2f',rRall)])
axis([0 1.1 0 1.1])
xlabel('Reg')
ylabel('POEM')

subplot(3,3,5)
plot(f_large(id4),Col_frac(id4),'.k','MarkerSize',15); hold on;
plot(x,x,'--r');hold on;
text(0.05,1,['r = ' sprintf('%2.2f',rRall4)])
axis([0 1.1 0 1.1])
xlabel('Reg')
ylabel('POEM')

subplot(3,3,6)
plot(f_large(id2),Col_frac(id2),'.k','MarkerSize',15); hold on;
plot(x,x,'--r');hold on;
text(0.05,1,['r = ' sprintf('%2.2f',rRall2)])
axis([0 1.1 0 1.1])
xlabel('Reg')
ylabel('POEM')
stamp(cfile)
print('-dpng',[ppath 'Clim_' harv '_Ecoreg_fracPD_catch_DvD_Reg_corr.png'])

