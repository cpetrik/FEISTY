% Visualize output of FEISTY forced with CESM
% Time series plots of Fish-MIP calcs

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'FishMIP_output_Preindust_' cfile '.mat']);
[ni,nj,nt] = size(tsb);
Ptsb=nanmean(reshape(tsb,ni*nj,nt),1);
Ptcb=nanmean(reshape(tcb,ni*nj,nt),1);
Pb10cm=nanmean(reshape(b10cm,ni*nj,nt),1);
Pb30cm=nanmean(reshape(b30cm,ni*nj,nt),1);
clear tsb tcb b10cm b30cm

load([fpath 'FishMIP_output_Historic_pristine_' cfile '.mat']);
[ni,nj,nt] = size(tsb);
Htsb=nanmean(reshape(tsb,ni*nj,nt),1);
Htcb=nanmean(reshape(tcb,ni*nj,nt),1);
Hb10cm=nanmean(reshape(b10cm,ni*nj,nt),1);
Hb30cm=nanmean(reshape(b30cm,ni*nj,nt),1);
clear tsb tcb b10cm b30cm

load([fpath 'FishMIP_output_Forecast_pristine_' cfile '.mat']);
[ni,nj,nt] = size(tsb);
Rtsb=nanmean(reshape(tsb,ni*nj,nt),1);
Rtcb=nanmean(reshape(tcb,ni*nj,nt),1);
Rb10cm=nanmean(reshape(b10cm,ni*nj,nt),1);
Rb30cm=nanmean(reshape(b30cm,ni*nj,nt),1);
clear tsb tcb b10cm b30cm

load([fpath 'FishMIP_output_NPP_cont_' cfile '.mat']);
[ni,nj,nt] = size(tsb);
Ntsb=nanmean(reshape(tsb,ni*nj,nt),1);
Ntcb=nanmean(reshape(tcb,ni*nj,nt),1);
Nb10cm=nanmean(reshape(b10cm,ni*nj,nt),1);
Nb30cm=nanmean(reshape(b30cm,ni*nj,nt),1);
clear tsb tcb b10cm b30cm

load([fpath 'FishMIP_output_Temp_cont_' cfile '.mat']);
[ni,nj,nt] = size(tsb);
Ttsb=nanmean(reshape(tsb,ni*nj,nt),1);
Ttcb=nanmean(reshape(tcb,ni*nj,nt),1);
Tb10cm=nanmean(reshape(b10cm,ni*nj,nt),1);
Tb30cm=nanmean(reshape(b30cm,ni*nj,nt),1);
clear tsb tcb b10cm b30cm

%% Plots in time
yrs = 1850:2100;
hyrs = 1850:2005;
fyrs = 2006:2100;

figure(1)
plot(yrs,log10(Ptsb),'k','Linewidth',2); hold on;
plot(hyrs,log10(Htsb),'--b','Linewidth',2); hold on;
plot(fyrs,log10(Rtsb),'--r','Linewidth',2); hold on;
plot(yrs,log10(Ntsb),':m','Linewidth',2); hold on;
plot(yrs,log10(Ttsb),':g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('log_1_0 Total System Biomass (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_tsb.png'])

%%
figure(2)
plot(yrs,log10(Ptcb),'k','Linewidth',2); hold on;
plot(hyrs,log10(Htcb),'--b','Linewidth',2); hold on;
plot(fyrs,log10(Rtcb),'--r','Linewidth',2); hold on;
plot(yrs,log10(Ntcb),':m','Linewidth',2); hold on;
plot(yrs,log10(Ttcb),':g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
title('log_1_0 Total Consumer Biomass (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_tcb.png'])

%%
figure(3)
plot(yrs,log10(Pb10cm),'k','Linewidth',2); hold on;
plot(hyrs,log10(Hb10cm),'b','Linewidth',2); hold on;
plot(fyrs,log10(Rb10cm),'r','Linewidth',2); hold on;
plot(yrs,log10(Nb10cm),'m','Linewidth',2); hold on;
plot(yrs,log10(Tb10cm),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
title('log_1_0 Consumers > 10cm Biomass (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_b10cm.png'])

%%
figure(4)
plot(yrs,log10(Pb30cm),'k','Linewidth',2); hold on;
plot(hyrs,log10(Hb30cm),'b','Linewidth',2); hold on;
plot(fyrs,log10(Rb30cm),'r','Linewidth',2); hold on;
plot(yrs,log10(Nb30cm),'m','Linewidth',2); hold on;
plot(yrs,log10(Tb30cm),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
title('log_1_0 Consumers > 30cm Biomass (gC m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_b30cm.png'])

%% Save
save([fpath 'FishMIP_output_timeseries_' cfile '.mat']);
