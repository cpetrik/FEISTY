% Compare daily interp zooplankton and detritus from
% CMIP6 COBALT output and ESM2M output
% to check unit conversion

clear all
close all

%% Hist
%Zoop and det
load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/Data_gfdl_hist_daily_2000.mat');
load('/Volumes/FEISTY/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2000.mat');

%%
%Is ESM off by 10? - NO
% test = ESM.Zm*10;

hzm1=nanmean(COBALT.Zm) + nanmean(COBALT.Zl);
hzm2=nanmean(ESM.Zm);
% hzm3=nanmean(test);

hd1=nanmean(COBALT.det);
hd2=nanmean(ESM.det);


figure(1)
plot(hzm1,'b'); hold on;
plot(hzm2,'r'); hold on;
% plot(hzm3,'k'); hold on;
title('Historic mesozoop')
legend('COB','ESM')%,'test')

%
figure(2)
plot(hd1,'k'); hold on;
plot(hd2,'color',[0.5 0.5 0.5]); hold on;
title('Historic det')
legend('COB','ESM')

% Comp Hist
hsperc = 100*(hzm1 - hzm2) ./ hzm1;

figure(3)
subplot(2,2,1)
histogram(hsperc)
title('Historic mesozoop percent diff')

%% RCP 85
%Zoop and det
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/ssp585/Data_gfdl_ssp585_daily_2050.mat');
load('/Volumes/FEISTY/POEM_JLD/rcp85/Data_rcp85_2050.mat');

fzm1=nanmean(COBALT.Zm) + nanmean(COBALT.Zl);
fzm2=nanmean(ESM.Zm);

fd1=nanmean(COBALT.det);
fd2=nanmean(ESM.det);

%
figure(4)
plot(fzm1,'b'); hold on;
plot(fzm2,'r'); hold on;
title('RCP 8.5 mesozoop')

%
figure(5)
plot(fd1,'k'); hold on;
plot(fd2,'color',[0.5 0.5 0.5]); hold on;
title('RCP 8.5 det')

% Comp
fsperc = 100*(fzm1 - fzm2) ./ fzm1;

figure(3)
subplot(2,2,2)
histogram(fsperc)
title('RCP 8.5 mesozoop percent diff')

%% Preindust
%Zoop and det
load('/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/preindust/Data_gfdl_spinup_daily_1860.mat');
load('/Volumes/FEISTY/POEM_JLD/pre_indust/Data_preindust_1860.mat');

pzm1=nanmean(COBALT.Zm) + nanmean(COBALT.Zl);
pzm2=nanmean(ESM.Zm);

pd1=nanmean(COBALT.det);
pd2=nanmean(ESM.det);

%
figure(6)
plot(pzm1,'b'); hold on;
plot(pzm2,'r'); hold on;
title('Preindustrial mesozoop')

%
figure(7)
plot(pd1,'k'); hold on;
plot(pd2,'color',[0.5 0.5 0.5]); hold on;
title('Preindustrial det')

% Comp
psperc = 100*(pzm1 - pzm2) ./ pzm1;

figure(3)
subplot(2,2,3)
histogram(psperc)
title('Preindustrial mesozoop percent diff')

