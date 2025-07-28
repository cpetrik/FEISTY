% Plot all scenarios together

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

esms = {'IPSL','UKESM','CESM2-WACCM','CESM2-WACCM'};
e2 = {'ipsl','ukesm','cesm2','cesm2'};

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/FishMIP/wg2300/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

fpath=['/Volumes/petrik-lab/Feisty/NC/WG2300/',cfile,'/IPSL/'];

%% Hist
load([fpath 'Means_IPSL_historic_pristine_' cfile '.mat'],...
    'mo',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean');

HistF_ts = sf_tmean + mf_tmean;
HistP_ts = sp_tmean + mp_tmean + lp_tmean;
HistD_ts = sd_tmean + md_tmean + ld_tmean;
HistB_ts = b_tmean;

HistA_ts = HistF_ts + HistP_ts + HistD_ts; 

Hist_time = mo;

clear mo sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean
clear lp_tmean ld_tmean b_tmean

%% SSP 126
load([fpath 'Means_IPSL_ssp126_pristine_' cfile '.mat'],...
    'mo',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean');

S126F_ts = sf_tmean + mf_tmean;
S126P_ts = sp_tmean + mp_tmean + lp_tmean;
S126D_ts = sd_tmean + md_tmean + ld_tmean;
S126B_ts = b_tmean;

S126A_ts = S126F_ts + S126P_ts + S126D_ts; 

S126_time = mo;

clear mo sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean
clear lp_tmean ld_tmean b_tmean);

%% SSP 585
load([fpath 'Means_IPSL_ssp585_pristine_' cfile '.mat'],...
    'mo',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean');

S585F_ts = sf_tmean + mf_tmean;
S585P_ts = sp_tmean + mp_tmean + lp_tmean;
S585D_ts = sd_tmean + md_tmean + ld_tmean;
S585B_ts = b_tmean;

S585A_ts = S585F_ts + S585P_ts + S585D_ts; 

S585_time = mo;

clear mo sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean
clear lp_tmean ld_tmean b_tmean);

%% SSP 534
load([fpath 'Means_IPSL_ssp534-over_pristine_' cfile '.mat'],...
    'mo',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean');

S534F_ts = sf_tmean + mf_tmean;
S534P_ts = sp_tmean + mp_tmean + lp_tmean;
S534D_ts = sd_tmean + md_tmean + ld_tmean;
S534B_ts = b_tmean;

S534A_ts = S534F_ts + S534P_ts + S534D_ts; 

S534_time = mo;

clear mo sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean
clear lp_tmean ld_tmean b_tmean);

%% Rolling means ~ annual
hist_yr = movmean(Hist_time,12);
ssp126_yr = movmean(S126_time,12);
ssp585_yr = movmean(S585_time,12);
ssp534_yr = movmean(S534_time,12);

hist_Tp = movmean(HistF_ts,12);
ssp126_Tp = movmean(S126F_ts,12);
ssp585_Tp = movmean(S585F_ts,12);
ssp534_Tp = movmean(S534F_ts,12);

hist_Tb = movmean(HistP_ts,12);
ssp126_Tb = movmean(S126P_ts,12);
ssp585_Tb = movmean(S585P_ts,12);
ssp534_Tb = movmean(S534P_ts,12);

hist_Zm = movmean(HistD_ts,12);
ssp126_Zm = movmean(S126D_ts,12);
ssp585_Zm = movmean(S585D_ts,12);
ssp534_Zm = movmean(S534D_ts,12);

hist_Det = movmean(HistA_ts,12);
ssp126_Det = movmean(S126A_ts,12);
ssp585_Det = movmean(S585A_ts,12);
ssp534_Det = movmean(S534A_ts,12);

%% biomass
figure(1)
subplot(2,2,1)
plot(hist_yr,hist_Tp,'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tp,'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tp,'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tp,'color',[0 0.75 0.5],'LineWidth',2);
title('IPSL Forage')
legend('Hist','126','585','534')
legend('location','southwest')

subplot(2,2,2)
plot(hist_yr,hist_Tb,'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tb,'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tb,'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tb,'color',[0 0.75 0.5],'LineWidth',2);
title('IPSL Lg Pel')

subplot(2,2,3)
plot(hist_yr,hist_Zm,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Zm,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Zm,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Zm,'color',[0 0.75 0.5],'LineWidth',1.5);
title('IPSL Dem')

subplot(2,2,4)
plot(hist_yr,hist_Det,'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Det,'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Det,'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Det,'color',[0 0.75 0.5],'LineWidth',1.5);
title('IPSL All fishes')
print('-dpng',[ppath 'IPSL_global_time_means_all_scenarios.png'])

%% biomass diffs
figure(2)
subplot(2,2,1)
plot(hist_yr,hist_Tp - hist_Tp(end),'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tp - hist_Tp(end),'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tp - hist_Tp(end),'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tp - hist_Tp(end),'color',[0 0.75 0.5],'LineWidth',2);
ylim([-0.5 0.2])
title('IPSL Forage')
legend('Hist','126','585','534')
legend('location','southwest')

subplot(2,2,2)
plot(hist_yr,hist_Tb - hist_Tb(end),'k','LineWidth',2); hold on
plot(ssp126_yr,ssp126_Tb - hist_Tb(end),'b','LineWidth',2); hold on
plot(ssp585_yr,ssp585_Tb - hist_Tb(end),'r','LineWidth',2); hold on
plot(ssp534_yr,ssp534_Tb - hist_Tb(end),'color',[0 0.75 0.5],'LineWidth',2);
ylim([-1 0.2])
title('IPSL Lg Pel')

subplot(2,2,3)
plot(hist_yr,hist_Zm - hist_Zm(end),'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Zm - hist_Zm(end),'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Zm - hist_Zm(end),'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Zm - hist_Zm(end),'color',[0 0.75 0.5],'LineWidth',1.5);
ylim([-0.5 0.2])
title('IPSL Dem')

subplot(2,2,4)
plot(hist_yr,hist_Det - hist_Det(end),'k','LineWidth',1.5); hold on
plot(ssp126_yr,ssp126_Det - hist_Det(end),'b','LineWidth',1.5); hold on
plot(ssp585_yr,ssp585_Det - hist_Det(end),'r','LineWidth',1.5); hold on
plot(ssp534_yr,ssp534_Det - hist_Det(end),'color',[0 0.75 0.5],'LineWidth',1.5);
ylim([-1.5 0.1])
title('IPSL All fishes')
print('-dpng',[ppath 'IPSL_global_time_diffs_all_scenarios.png'])


%% biomass - non rolling means b/c SSP 534 looks off
figure(3)
subplot(2,2,1)
plot(Hist_time,HistF_ts,'k','LineWidth',2); hold on
plot(S126_time,S126F_ts,'b','LineWidth',2); hold on
plot(S585_time,S585F_ts,'r','LineWidth',2); hold on
plot(S534_time,S534F_ts,'color',[0 0.75 0.5],'LineWidth',2);
title('IPSL Forage')
legend('Hist','126','585','534')
legend('location','southwest')

subplot(2,2,2)
plot(Hist_time,HistP_ts,'k','LineWidth',2); hold on
plot(S126_time,S126P_ts,'b','LineWidth',2); hold on
plot(S585_time,S585P_ts,'r','LineWidth',2); hold on
plot(S534_time,S534P_ts,'color',[0 0.75 0.5],'LineWidth',2);
title('IPSL Lg Pel')

subplot(2,2,3)
plot(Hist_time,HistD_ts,'k','LineWidth',1.5); hold on
plot(S126_time,S126D_ts,'b','LineWidth',1.5); hold on
plot(S585_time,S585D_ts,'r','LineWidth',1.5); hold on
plot(S534_time,S534D_ts,'color',[0 0.75 0.5],'LineWidth',1.5);
title('IPSL Dem')

subplot(2,2,4)
plot(Hist_time,HistA_ts,'k','LineWidth',1.5); hold on
plot(S126_time,S126A_ts,'b','LineWidth',1.5); hold on
plot(S585_time,S585A_ts,'r','LineWidth',1.5); hold on
plot(S534_time,S534A_ts,'color',[0 0.75 0.5],'LineWidth',1.5);
title('IPSL All fishes')
%print('-dpng',[ppath 'IPSL_global_time_means_all_scenarios.png'])

%% test 
load(['/Volumes/petrik-lab/Feisty/NC/WG2300/',cfile,'/IPSL/Last_mo_IPSL_ssp585_pristine_',cfile,'.mat']);

S534F1 = mean((Sml_f.bio + Med_f.bio),'omitnan') %0.574
S534F_ts(1) %similar                             %0.579
S534P1 = mean((Sml_p.bio + Med_p.bio + Lrg_p.bio),'omitnan') %1.297
S534P_ts(1) %way different                                   %0.4853
S534D1 = mean((Sml_d.bio + Med_d.bio + Lrg_d.bio),'omitnan') %1.774
S534D_ts(1)                                                  %1.522
S534B1 = mean(BENT.bio,'omitnan') %0.6347
S534B_ts(1)                       %1.091

S534A1 = S534F1 + S534P1 + S534D1 % 3.6413
S534A_ts(1)                       % 2.5858




