%% Check that daily interp files got made correctly

clear all
close all

% Preindust
load(['/Volumes/GFDL/Fish-MIP/CESM/PreIndust/Data_cesm_pi_2005.mat']);
PI1 = CESM;
clear CESM
load(['/Volumes/GFDL/Fish-MIP/CESM/PreIndust/Data_cesm_pi_2100.mat']);
PI2 = CESM;
clear CESM
% Hist
load(['/Volumes/GFDL/Fish-MIP/CESM/Hist/Data_cesm_hist_2005.mat']);
HI = CESM;
clear CESM
% RCP 85
load(['/Volumes/GFDL/Fish-MIP/CESM/RCP85/Data_cesm_rcp85_2100.mat']);
RC = CESM;
clear CESM
% Temp cont
load(['/Volumes/GFDL/Fish-MIP/CESM/Temp_cont/Data_cesm_temp_cont_2005.mat']);
TC1 = CESM;
clear CESM
load(['/Volumes/GFDL/Fish-MIP/CESM/Temp_cont/Data_cesm_temp_cont_2100.mat']);
TC2 = CESM;
clear CESM
% NPP cont
load(['/Volumes/GFDL/Fish-MIP/CESM/NPP_cont/Data_cesm_npp_cont_2005.mat']);
NP1 = CESM;
clear CESM
load(['/Volumes/GFDL/Fish-MIP/CESM/NPP_cont/Data_cesm_npp_cont_2100.mat']);
NP2 = CESM;
clear CESM

%%
zm_pi1=nanmean(PI1.Zm);
zl_pi1=nanmean(PI1.Zl);
zm_pi2=nanmean(PI2.Zm);
zl_pi2=nanmean(PI2.Zl);

zm_hi=nanmean(HI.Zm);
zl_hi=nanmean(HI.Zl);

zm_rc=nanmean(RC.Zm);
zl_rc=nanmean(RC.Zl);

zm_tp1=nanmean(TC1.Zm);
zl_tp1=nanmean(TC1.Zl);
zm_tp2=nanmean(TC2.Zm);
zl_tp2=nanmean(TC2.Zl);

zm_np1=nanmean(NP1.Zm);
zl_np1=nanmean(NP1.Zl);
zm_np2=nanmean(NP2.Zm);
zl_np2=nanmean(NP2.Zl);

%% det
d_pi1=nanmean(PI1.det);
d_pi2=nanmean(PI2.det);

d_hi=nanmean(HI.det);

d_rc=nanmean(RC.det);

d_tp1=nanmean(TC1.det);
d_tp2=nanmean(TC2.det);

d_np1=nanmean(NP1.det);
d_np2=nanmean(NP2.det);

%% temp
tp_pi1=nanmean(PI1.Tp);
tb_pi1=nanmean(PI1.Tb);
tp_pi2=nanmean(PI2.Tp);
tb_pi2=nanmean(PI2.Tb);

tp_hi=nanmean(HI.Tp);
tb_hi=nanmean(HI.Tb);

tp_rc=nanmean(RC.Tp);
tb_rc=nanmean(RC.Tb);

tp_tp1=nanmean(TC1.Tp);
tb_tp1=nanmean(TC1.Tb);
tp_tp2=nanmean(TC2.Tp);
tb_tp2=nanmean(TC2.Tb);

tp_np1=nanmean(NP1.Tp);
tb_np1=nanmean(NP1.Tb);
tp_np2=nanmean(NP2.Tp);
tb_np2=nanmean(NP2.Tb);

%% 2005
%Det
figure(1)
plot(d_pi1,'k'); hold on;
plot(d_hi,'b'); hold on;
plot(d_tp1,'--c'); hold on;
plot(d_np1,'--m'); hold on;
title('2005 det')
legend({'pi','hist','tempc','nppc'})

%M Zoop
figure(2)
plot(zm_pi1,'k'); hold on;
plot(zm_hi,'b'); hold on;
plot(zm_tp1,'--c'); hold on;
plot(zm_np1,'--m'); hold on;
title('2005 MZ')
legend({'pi','hist','tempc','nppc'})

%L Zoop
figure(3)
plot(zl_pi1,'k'); hold on;
plot(zl_hi,'b'); hold on;
plot(zl_tp1,'--c'); hold on;
plot(zl_np1,'--m'); hold on;
title('2005 LZ')
legend({'pi','hist','tempc','nppc'})

%Temp P
figure(4)
plot(tp_pi1,'k'); hold on;
plot(tp_hi,'b'); hold on;
plot(tp_tp1,'--c'); hold on;
plot(tp_np1,'--m'); hold on;
title('2005 Tp')
legend({'pi','hist','tempc','nppc'})

%Temp B
figure(5)
plot(tb_pi1,'k'); hold on;
plot(tb_hi,'b'); hold on;
plot(tb_tp1,'--c'); hold on;
plot(tb_np1,'--m'); hold on;
title('2005 Tb')
legend({'pi','hist','tempc','nppc'})

%% 2100
%Det
figure(6)
plot(d_pi2,'k'); hold on;
plot(d_rc,'r'); hold on;
plot(d_tp2,'--c'); hold on;
plot(d_np2,'--m'); hold on;
title('2100 det')
legend({'pi','rcp','tempc','nppc'})

%M Zoop
figure(7)
plot(zm_pi2,'k'); hold on;
plot(zm_rc,'r'); hold on;
plot(zm_tp2,'--c'); hold on;
plot(zm_np2,'--m'); hold on;
title('2100 MZ')
legend({'pi','rcp','tempc','nppc'})

%L Zoop
figure(8)
plot(zl_pi2,'k'); hold on;
plot(zl_rc,'r'); hold on;
plot(zl_tp2,'--c'); hold on;
plot(zl_np2,'--m'); hold on;
title('2100 LZ')
legend({'pi','rcp','tempc','nppc'})

%Temp P
figure(9)
plot(tp_pi2,'k'); hold on;
plot(tp_rc,'r'); hold on;
plot(tp_tp2,'--c'); hold on;
plot(tp_np2,'--m'); hold on;
title('2100 Tp')
legend({'pi','rcp','tempc','nppc'})

%Temp B
figure(10)
plot(tb_pi2,'k'); hold on;
plot(tb_rc,'r'); hold on;
plot(tb_tp2,'--c'); hold on;
plot(tb_np2,'--m'); hold on;
title('2100 Tb')
legend({'pi','rcp','tempc','nppc'})



