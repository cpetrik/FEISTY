% Visualize output of FEISTY forced with CESM
% Time series plots of Fish-MIP calcs

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'NoDc_enc70-b200_m4-b175-k086_c20-b250_noD_J100_A050_Sm025_nmort1_BE00_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];
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

%%
load([fpath 'Ts_Means_Preindust_' cfile '.mat']);
load([fpath 'Ts_Means_Historic_' cfile '.mat']);
load([fpath 'Ts_Means_Forecast_' cfile '.mat']);
load([fpath 'Ts_Means_NPP_cont_' cfile '.mat']);
load([fpath 'Ts_Means_Temp_cont_' cfile '.mat']);

PA = PF + PP;
HA = HF + HP;
RA = RF + RP;
NA = NF + NP;
TA = TF + TP;

%% Relative to 1860-1870 mean
yid = find(Py<1870 & Py>=1860);

rPF = PF - mean(PF(yid));
rPP = PP - mean(PP(yid));
rPA = PA - mean(PA(yid));

rHF = HF - mean(HF(yid));
rHP = HP - mean(HP(yid));
rHA = HA - mean(HA(yid));

rNF = NF - mean(NF(yid));
rNP = NP - mean(NP(yid));
rNA = NA - mean(NA(yid));

rTF = TF - mean(TF(yid));
rTP = TP - mean(TP(yid));
rTA = TA - mean(TA(yid));

rRF = RF - mean(HF(yid));
rRP = RP - mean(HP(yid));
rRA = RA - mean(HA(yid));

%% Plots in time
yrs = 1850:2100;
hyrs = 1850:2005;
fyrs = 2006:2100;

%%
figure(10)
plot(Py,log10(PF),'k','Linewidth',2); hold on;
plot(Hy,log10(HF),'b','Linewidth',2); hold on;
plot(Ry,log10(RF),'r','Linewidth',2); hold on;
plot(Ny,log10(NF),'m','Linewidth',2); hold on;
plot(Ty,log10(TF),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('log_1_0 F Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_F.png'])

figure(11)
plot(Py,log10(PP),'k','Linewidth',2); hold on;
plot(Hy,log10(HP),'b','Linewidth',2); hold on;
plot(Ry,log10(RP),'r','Linewidth',2); hold on;
plot(Ny,log10(NP),'m','Linewidth',2); hold on;
plot(Ty,log10(TP),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('log_1_0 P Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_P.png'])

figure(12)
plot(Py,log10(PA),'k','Linewidth',2); hold on;
plot(Hy,log10(HA),'b','Linewidth',2); hold on;
plot(Ry,log10(RA),'r','Linewidth',2); hold on;
plot(Ny,log10(NA),'m','Linewidth',2); hold on;
plot(Ty,log10(TA),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('log_1_0 Fish Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_All.png'])


figure(13)
plot(Py,(rPF),'k','Linewidth',2); hold on;
plot(Hy,(rHF),'b','Linewidth',2); hold on;
plot(Ry,(rRF),'r','Linewidth',2); hold on;
plot(Ny,(rNF),'m','Linewidth',2); hold on;
plot(Ty,(rTF),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('F Biomass (g m^-^2) relative to 1860-1870')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_F.png'])

figure(14)
plot(Py,(rPP),'k','Linewidth',2); hold on;
plot(Hy,(rHP),'b','Linewidth',2); hold on;
plot(Ry,(rRP),'r','Linewidth',2); hold on;
plot(Ny,(rNP),'m','Linewidth',2); hold on;
plot(Ty,(rTP),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('P Biomass (g m^-^2) relative to 1860-1870')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_P.png'])

figure(15)
plot(Py,(rPA),'k','Linewidth',2); hold on;
plot(Hy,(rHA),'b','Linewidth',2); hold on;
plot(Ry,(rRA),'r','Linewidth',2); hold on;
plot(Ny,(rNA),'m','Linewidth',2); hold on;
plot(Ty,(rTA),'g','Linewidth',2); hold on;
xlim([yrs(1) yrs(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
xlabel('Year')
title('Fish Biomass (g m^-^2) relative to 1860-1870')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_All.png'])


%%
figure(1)
plot(yrs,log10(Ptsb),'k','Linewidth',2); hold on;
plot(hyrs,log10(Htsb),'b','Linewidth',2); hold on;
plot(fyrs,log10(Rtsb),'r','Linewidth',2); hold on;
plot(yrs,log10(Ntsb),'m','Linewidth',2); hold on;
plot(yrs,log10(Ttsb),'g','Linewidth',2); hold on;
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
plot(hyrs,log10(Htcb),'b','Linewidth',2); hold on;
plot(fyrs,log10(Rtcb),'r','Linewidth',2); hold on;
plot(yrs,log10(Ntcb),'m','Linewidth',2); hold on;
plot(yrs,log10(Ttcb),'g','Linewidth',2); hold on;
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
