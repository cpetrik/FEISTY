% Visualize output of FEISTY forced with CESM
% Preindustrial 1800-2100 initialized with spinup biomass
% Time series plots and maps
% 2090-2100_v_1860-1870
% Change biomass to relative to 1860-1870 mean

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

%% Relative to 1860-1870 mean
yid = find(Py<1870 & Py>=1860);

rPF = PF - mean(PF(yid));
rPP = PP - mean(PP(yid));
rPD = PD - mean(PD(yid));
rPB = PB - mean(PB(yid));

rHF = HF - mean(HF(yid));
rHP = HP - mean(HP(yid));
rHD = HD - mean(HD(yid));
rHB = HB - mean(HB(yid));

rNF = NF - mean(NF(yid));
rNP = NP - mean(NP(yid));
rND = ND - mean(ND(yid));
rNB = NB - mean(NB(yid));

rTF = TF - mean(TF(yid));
rTP = TP - mean(TP(yid));
rTD = TD - mean(TD(yid));
rTB = TB - mean(TB(yid));

rRF = RF - mean(HF(yid));
rRP = RP - mean(HP(yid));
rRD = RD - mean(HD(yid));
rRB = RB - mean(HB(yid));

%% Plots in time
figure(1)
%plot(Py,(PB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(Py,(rPF),'k','Linewidth',2); hold on;
plot(Py,(rPP),'k','Linewidth',2); hold on;
plot(Py,(rPD),'k','Linewidth',2); hold on;

plot(Hy,(rHF),'b','Linewidth',2); hold on;
plot(Hy,(rHP),'b','Linewidth',2); hold on;
plot(Hy,(rHD),'b','Linewidth',2); hold on;

plot(Ry,(rRF),'r','Linewidth',2); hold on;
plot(Ry,(rRP),'r','Linewidth',2); hold on;
plot(Ry,(rRD),'r','Linewidth',2); hold on;

plot(Ny,(rNF),'m','Linewidth',2); hold on;
plot(Ny,(rNP),'m','Linewidth',2); hold on;
plot(Ny,(rND),'m','Linewidth',2); hold on;

plot(Ty,(rTF),'g','Linewidth',2); hold on;
plot(Ty,(rTP),'g','Linewidth',2); hold on;
plot(Ty,(rTD),'g','Linewidth',2); hold on;
%legend('F','P','D')
%legend('location','northeast')
xlim([Py(1) Py(end)])
%ylim([-1 1])
xlabel('Year')
ylabel('Biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_all_types.png'])

%%
figure(2)
plot(Py,(rPF),'k','Linewidth',2); hold on;
plot(Hy,(rHF),'b','Linewidth',2); hold on;
plot(Ry,(rRF),'r','Linewidth',2); hold on;
plot(Ny,(rNF),'m','Linewidth',2); hold on;
plot(Ty,(rTF),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPPcon','Tempcon')
legend('location','northwest')
title(['Forage fish ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_F.png'])

figure(3)
plot(Py,(rPD),'k','Linewidth',2); hold on;
plot(Hy,(rHD),'b','Linewidth',2); hold on;
plot(Ry,(rRD),'r','Linewidth',2); hold on;
plot(Ny,(rND),'m','Linewidth',2); hold on;
plot(Ty,(rTD),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPPcon','Tempcon')
legend('location','northwest')
title(['Demersal fish ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_D.png'])

figure(4)
plot(Py,(rPB),'k','Linewidth',2); hold on;
plot(Hy,(rHB),'b','Linewidth',2); hold on;
plot(Ry,(rRB),'r','Linewidth',2); hold on;
plot(Ny,(rNB),'m','Linewidth',2); hold on;
plot(Ty,(rTB),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPPcon','Tempcon')
legend('location','northwest')
xlabel('Year')
title(['Benthic invertebrates ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_B.png'])
%%
figure(5)
plot(Py,(rPP),'k','Linewidth',2); hold on;
plot(Hy,(rHP),'b','Linewidth',2); hold on;
plot(Ry,(rRP),'r','Linewidth',2); hold on;
plot(Ny,(rNP),'m','Linewidth',2); hold on;
plot(Ty,(rTP),'g','Linewidth',2); hold on;
xlim([Py(1) Py(end)])
legend('Pre','Hist','RCP','NPPcon','Tempcon')
legend('location','southwest')
title(['Large pelagic fish ' harv ' biomass (g m^-^2)'])
stamp(cfile)
print('-dpng',[ppath 'Timeseries_rel_',harv,'_P.png'])


