clear all
close all

%% Hist
%Zoop and det
load(['/Volumes/GFDL/Fish-MIP/CESM/Hist/Data_cesm_hist_2000.mat']);
load(['/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2000.mat']);

hzm1=nanmean(COBALT.Zm);
hzl1=nanmean(COBALT.Zl);
hzm2=nanmean(CESM.Zm);
hzl2=nanmean(CESM.Zl);

hd1=nanmean(COBALT.det);
hd2=nanmean(CESM.det);

%
figure(1)
plot(hzm1,'b'); hold on;
plot(hzl1,'r'); hold on;
plot(1e3*hzm2,'color',[0 0 0.5]); hold on;
plot(1e3*hzl2,'color',[0.5 0 0]); hold on;
title('Historic')

%
figure(2)
plot(hd1,'k'); hold on;
plot(hd2,'color',[0.5 0.5 0.5]); hold on;
title('Historic')

% Comp Hist
hsperc = 100*(hzm1 - 1e3*hzm2) ./ hzm1;
hlperc = 100*(hzl1 - 1e3*hzl2) ./ hzl1;

figure(9)
subplot(1,2,1)
hist(hsperc)
subplot(1,2,2)
hist(hlperc)
title('Historic')

%% RCP 85
%Zoop and det
load(['/Volumes/GFDL/Fish-MIP/CESM/RCP85/Data_cesm_rcp85_2050.mat']);
load(['/Volumes/GFDL/POEM_JLD/rcp85/Data_rcp85_2050.mat']);

fzm1=nanmean(COBALT.Zm);
fzl1=nanmean(COBALT.Zl);
fzm2=nanmean(CESM.Zm);
fzl2=nanmean(CESM.Zl);

fd1=nanmean(COBALT.det);
fd2=nanmean(CESM.det);

%
figure(3)
plot(fzm1,'b'); hold on;
plot(fzl1,'r'); hold on;
plot(1e3*fzm2,'color',[0 0 0.5]); hold on;
plot(1e3*fzl2,'color',[0.5 0 0]); hold on;
title('RCP 8.5')

%
figure(4)
plot(fd1,'k'); hold on;
plot(fd2,'color',[0.5 0.5 0.5]); hold on;
title('RCP 8.5')

% Comp
fsperc = 100*(fzm1 - 1e3*fzm2) ./ fzm1;
flperc = 100*(fzl1 - 1e3*fzl2) ./ fzl1;

figure(8)
subplot(1,2,1)
hist(fsperc)
subplot(1,2,2)
hist(flperc)
title('RCP 8.5')

%% Preindust
%Zoop and det
load(['/Volumes/GFDL/Fish-MIP/CESM/PreIndust/Data_cesm_pi_1860.mat']);
load(['/Volumes/GFDL/POEM_JLD/pre_indust/Data_preindust_1860.mat']);

pzm1=nanmean(COBALT.Zm);
pzl1=nanmean(COBALT.Zl);
pzm2=nanmean(CESM.Zm);
pzl2=nanmean(CESM.Zl);

pd1=nanmean(COBALT.det);
pd2=nanmean(CESM.det);

%
figure(5)
plot(pzm1,'b'); hold on;
plot(pzl1,'r'); hold on;
plot(1e3*pzm2,'color',[0 0 0.5]); hold on;
plot(1e3*pzl2,'color',[0.5 0 0]); hold on;
title('Preindustrial')

%
figure(6)
plot(pd1,'k'); hold on;
plot(pd2,'color',[0.5 0.5 0.5]); hold on;
title('Preindustrial')

% Comp
psperc = 100*(pzm1 - 1e3*pzm2) ./ pzm1;
plperc = 100*(pzl1 - 1e3*pzl2) ./ pzl1;

figure(7)
subplot(1,2,1)
hist(psperc)
subplot(1,2,2)
hist(plperc)
title('Preindustrial')

