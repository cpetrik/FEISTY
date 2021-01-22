% Check that daily interp files got made correctly
% And that units are correct

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% CORE
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_ocean_cobalt_daily_2005.mat');
CI = COBALT;
clear COBALT

%% Hist
load('/Volumes/FEISTY/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2005.mat');
HI = COBALT;
clear COBALT

%%
zm_ci=nanmean(CI.Zm);
zl_ci=nanmean(CI.Zl);
d_ci=nanmean(CI.det);
tp_ci=nanmean(CI.Tp);
tb_ci=nanmean(CI.Tb);
dzm_ci=nanmean(CI.dZm);
dzl_ci=nanmean(CI.dZl);

zm_hi=nanmean(HI.Zm);
zl_hi=nanmean(HI.Zl);
d_hi=nanmean(HI.det);
tp_hi=nanmean(HI.Tp);
tb_hi=nanmean(HI.Tb);
dzm_hi=nanmean(HI.dZm);
dzl_hi=nanmean(HI.dZl);

%% 2005
%Det
figure(1)
subplot(2,2,1)
plot(d_ci,'k'); hold on;
plot(d_hi,'b'); hold on;
title('2005 det')
legend({'core','hist'})

%Temp P
subplot(2,2,3)
plot(tp_ci,'k'); hold on;
plot(tp_hi,'b'); hold on;
title('2005 Tp')
%legend({'core','hist'})

%Temp B
subplot(2,2,4)
plot(tb_ci,'k'); hold on;
plot(tb_hi,'b'); hold on;
title('2005 Tb')
%legend({'core','hist'})
stamp('')
%print('-dpng',[pp 'GFDL_daily_ts_2005.png'])

%% 
%MZ
figure(2)
subplot(2,2,1)
plot(zm_ci,'k'); hold on;
plot(zm_hi,'b'); hold on;
title('2005 MZ')
legend({'core','hist'})

subplot(2,2,2)
plot(dzm_ci,'k'); hold on;
plot(dzm_hi,'b'); hold on;
title('2005 hplMZ')

%Temp P
subplot(2,2,3)
plot(zl_ci,'k'); hold on;
plot(zl_hi,'b'); hold on;
title('2005 LZ')
%legend({'core','hist'})

%Temp B
subplot(2,2,4)
plot(dzl_ci,'k'); hold on;
plot(dzl_hi,'b'); hold on;
title('2005 hplLZ')
%legend({'core','hist'})
stamp('')
%print('-dpng',[pp 'GFDL_daily_ts_2005.png'])

%%
%MZ
figure(3)
subplot(2,2,1)
plot(log10(zm_ci),'k'); hold on;
plot(log10(zm_hi),'b'); hold on;
title('2005 log10 MZ')
legend({'core','hist'})


%Temp P
subplot(2,2,3)
plot(log10(zl_ci),'k'); hold on;
plot(log10(zl_hi),'b'); hold on;
title('2005 log10 LZ')
%legend({'core','hist'})

%% Core zoo too low by 10^3 - Fixed now 
mean(HI.Zl(:))/mean(CI.Zl(:))
mean(HI.Zm(:))/mean(CI.Zm(:))


