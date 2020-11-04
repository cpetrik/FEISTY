%% Check that daily interp files got made correctly

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';

% Preindust
load(['/Volumes/MIP/Fish-MIP/CMIP6/GFDL/preindust/Data_gfdl_pi_daily_2005.mat']);
PI1 = ESM;
clear ESM
load(['/Volumes/MIP/Fish-MIP/CMIP6/GFDL/preindust/Data_gfdl_pi_daily_2099.mat']);
PI2 = ESM;
clear ESM
% Hist
load(['/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/Data_gfdl_hist_daily_2005.mat']);
HI = ESM;
clear ESM
% SSP 126
load(['/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp126/Data_gfdl_ssp126_daily_2099.mat']);
TC = ESM;
clear ESM
% RCP 585
load(['/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/Data_gfdl_ssp585_daily_2099.mat']);
RC = ESM;
clear ESM

%%
zm_pi1=nanmean(PI1.Zm);
zm_pi2=nanmean(PI2.Zm);

zm_hi=nanmean(HI.Zm);

zm_rc=nanmean(RC.Zm);

zm_tp=nanmean(TC.Zm);

%% det
d_pi1=nanmean(PI1.det);
d_pi2=nanmean(PI2.det);

d_hi=nanmean(HI.det);

d_rc=nanmean(RC.det);

d_tp=nanmean(TC.det);

%% temp
tp_pi1=nanmean(PI1.Tp);
tb_pi1=nanmean(PI1.Tb);
tp_pi2=nanmean(PI2.Tp);
tb_pi2=nanmean(PI2.Tb);

tp_hi=nanmean(HI.Tp);
tb_hi=nanmean(HI.Tb);

tp_rc=nanmean(RC.Tp);
tb_rc=nanmean(RC.Tb);

tp_tp=nanmean(TC.Tp);
tb_tp=nanmean(TC.Tb);

%% 2005
%Det
figure(1)
subplot(2,2,1)
plot(d_pi1,'k'); hold on;
plot(d_hi,'b'); hold on;
title('GFDL 2005 det')
legend({'pi','hist'})

%M Zoop
subplot(2,2,2)
plot(zm_pi1,'k'); hold on;
plot(zm_hi,'b'); hold on;
title('GFDL 2005 MZ')
%legend({'pi','hist'})

%Temp P
subplot(2,2,3)
plot(tp_pi1,'k'); hold on;
plot(tp_hi,'b'); hold on;
title('GFDL 2005 Tp')
%legend({'pi','hist'})

%Temp B
subplot(2,2,4)
plot(tb_pi1,'k'); hold on;
plot(tb_hi,'b'); hold on;
title('GFDL 2005 Tb')
%legend({'pi','hist'})
stamp('')
print('-dpng',[pp 'GFDL_daily_ts_2005.png'])

%% 2100
%Det
figure(2)
subplot(2,2,1)
plot(d_pi2,'k'); hold on;
plot(d_rc,'r'); hold on;
plot(d_tp,'m'); hold on;
title('GFDL 2100 det')
legend({'pi','585','126'})

%M Zoop
subplot(2,2,2)
plot(zm_pi2,'k'); hold on;
plot(zm_rc,'r'); hold on;
plot(zm_tp,'m'); hold on;
title('GFDL 2100 MZ')
%legend({'pi','585','126'})

%Temp P
subplot(2,2,3)
plot(tp_pi2,'k'); hold on;
plot(tp_rc,'r'); hold on;
plot(tp_tp,'m'); hold on;
title('GFDL 2100 Tp')
%legend({'pi','585','126'})

%Temp B
subplot(2,2,4)
plot(tb_pi2,'k'); hold on;
plot(tb_rc,'r'); hold on;
plot(tb_tp,'m'); hold on;
title('GFDL 2100 Tb')
%legend({'pi','585','126'})
stamp('')
print('-dpng',[pp 'GFDL_daily_ts_2100.png'])



