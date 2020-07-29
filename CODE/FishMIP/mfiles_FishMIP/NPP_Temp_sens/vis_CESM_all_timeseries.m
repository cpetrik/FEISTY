% Visualize output of FEISTY forced with CESM
% Preindustrial 1800-2100 initialized with spinup biomass
% Time series plots and maps

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Ts_Means_Preindust_' cfile '.mat']);
load([fpath 'Ts_Means_Historic_' cfile '.mat']);
load([fpath 'Ts_Means_Forecast_' cfile '.mat']);
load([fpath 'Ts_Means_NPP_cont_' cfile '.mat']);
load([fpath 'Ts_Means_Temp_cont_' cfile '.mat']);


%% Plots in time
figure(1)
%plot(Py,log10(PB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(Py,log10(PF),'k','Linewidth',2); hold on;
plot(Py,log10(PP),'k','Linewidth',2); hold on;
plot(Py,log10(PD),'k','Linewidth',2); hold on;

plot(Hy,log10(HF),'b','Linewidth',2); hold on;
plot(Hy,log10(HP),'b','Linewidth',2); hold on;
plot(Hy,log10(HD),'b','Linewidth',2); hold on;

plot(Ry,log10(RF),'r','Linewidth',2); hold on;
plot(Ry,log10(RP),'r','Linewidth',2); hold on;
plot(Ry,log10(RD),'r','Linewidth',2); hold on;

plot(Ny,log10(NF),'m','Linewidth',2); hold on;
plot(Ny,log10(NP),'m','Linewidth',2); hold on;
plot(Ny,log10(ND),'m','Linewidth',2); hold on;

plot(Ty,log10(TF),'g','Linewidth',2); hold on;
plot(Ty,log10(TP),'g','Linewidth',2); hold on;
plot(Ty,log10(TD),'g','Linewidth',2); hold on;
%legend('F','P','D')
%legend('location','northeast')
xlim([Py(1) Py(end)])
ylim([-1 1])
xlabel('Year')
ylabel('log_1_0 Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_all_types.png'])

%%
figure(2)
plot(Py,log10(PF),'k','Linewidth',2); hold on;
plot(Hy,log10(HF),'b','Linewidth',2); hold on;
plot(Ry,log10(RF),'r','Linewidth',2); hold on;
plot(Ny,log10(NF),'m','Linewidth',2); hold on;
plot(Ty,log10(TF),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','northwest')
title(['Forage fish ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_F.png'])

figure(3)
plot(Py,log10(PD),'k','Linewidth',2); hold on;
plot(Hy,log10(HD),'b','Linewidth',2); hold on;
plot(Ry,log10(RD),'r','Linewidth',2); hold on;
plot(Ny,log10(ND),'m','Linewidth',2); hold on;
plot(Ty,log10(TD),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','northwest')
title(['Demersal fish ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_D.png'])

figure(4)
plot(Py,log10(PB),'k','Linewidth',2); hold on;
plot(Hy,log10(HB),'b','Linewidth',2); hold on;
plot(Ry,log10(RB),'r','Linewidth',2); hold on;
plot(Ny,log10(NB),'m','Linewidth',2); hold on;
plot(Ty,log10(TB),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','northwest')
xlabel('Year')
title(['Benthic invertebrates ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_B.png'])
%%
figure(5)
plot(Py,log10(PP),'k','Linewidth',2); hold on;
plot(Hy,log10(HP),'b','Linewidth',2); hold on;
plot(Ry,log10(RP),'r','Linewidth',2); hold on;
plot(Ny,log10(NP),'m','Linewidth',2); hold on;
plot(Ty,log10(TP),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPP','Temp')
legend('location','southwest')
title(['Large pelagic fish ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_',harv,'_P.png'])


