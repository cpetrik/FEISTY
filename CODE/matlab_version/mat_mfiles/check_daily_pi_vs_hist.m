% Check that daily interp files got made correctly
% And that units are correct

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% CORE
load('/Volumes/MIP/GCM_DATA/ESM4_PI/Data_ocean_cobalt_daily_1400.mat');
CI = COBALT;
clear COBALT

%% Hist
load('/Volumes/FEISTY/POEM_JLD/esm2m_hist/Data_ESM2Mhist_1861.mat');
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

%% 
%Det
figure(1)
subplot(2,2,1)
plot(d_ci,'k'); hold on;
plot(d_hi,'b'); hold on;
title('det')
legend({'pi','hist'})

%Temp P
subplot(2,2,3)
plot(tp_ci,'k'); hold on;
plot(tp_hi,'b'); hold on;
title('Tp')
%legend({'pi','hist'})

%Temp B
subplot(2,2,4)
plot(tb_ci,'k'); hold on;
plot(tb_hi,'b'); hold on;
title('Tb')
%legend({'pi','hist'})
stamp('')
%print('-dpng',[pp 'GFDL_daily_ts_2005.png'])

%% 
%MZ
figure(2)
subplot(2,2,1)
plot(zm_ci,'k'); hold on;
plot(zm_hi,'b'); hold on;
title('MZ')
legend({'pi','hist'})

subplot(2,2,2)
plot(dzm_ci,'k'); hold on;
plot(dzm_hi,'b'); hold on;
title('hplMZ')

%Temp P
subplot(2,2,3)
plot(zl_ci,'k'); hold on;
plot(zl_hi,'b'); hold on;
title('LZ')
%legend({'pi','hist'})

%Temp B
subplot(2,2,4)
plot(dzl_ci,'k'); hold on;
plot(dzl_hi,'b'); hold on;
title('hplLZ')
%legend({'pi','hist'})
stamp('')
%print('-dpng',[pp 'GFDL_daily_ts_2005.png'])

%%
%MZ
figure(3)
subplot(2,2,1)
plot(log10(zm_ci),'k'); hold on;
plot(log10(zm_hi),'b'); hold on;
title('2005 log10 MZ')
legend({'pi','hist'})


%Temp P
subplot(2,2,3)
plot(log10(zl_ci),'k'); hold on;
plot(log10(zl_hi),'b'); hold on;
title('2005 log10 LZ')
%legend({'pi','hist'})



