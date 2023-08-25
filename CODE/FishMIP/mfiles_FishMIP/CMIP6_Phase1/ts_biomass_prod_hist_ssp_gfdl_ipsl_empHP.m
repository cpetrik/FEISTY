% Visualize time series output of FEISTY forced by CMIP6
% tcb relative to 1990-2000
% compared to tc production for Vianney's paper

clear 
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/FishMIP6/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% gfdl
gpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_CMIP6/' cfile '/'];
load([gpath 'Means_Hist_empHP_2000-2010_' cfile '.mat'],...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod');

%all consumers
GHtcb = sf_tmean+sp_tmean+sd_tmean+mf_tmean+mp_tmean+md_tmean+lp_tmean+ld_tmean;
GHtcp = sf_tprod+sp_tprod+sd_tprod+mf_tprod+mp_tprod+md_tprod+lp_tprod+ld_tprod;

%%
load([gpath 'Means_SSP585_empHP_2090-2100_' cfile '.mat'],...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod');

GStcb = sf_tmean+sp_tmean+sd_tmean+mf_tmean+mp_tmean+md_tmean+lp_tmean+ld_tmean;
GStcp = sf_tprod+sp_tprod+sd_tprod+mf_tprod+mp_tprod+md_tprod+lp_tprod+ld_tprod;

%% ipsl
ipath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/IPSL_CMIP6/' cfile '/'];
load([ipath 'Means_Hist_empHP_2000-2010_' cfile '.mat'],...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod');

IHtcb = sf_tmean+sp_tmean+sd_tmean+mf_tmean+mp_tmean+md_tmean+lp_tmean+ld_tmean;
IHtcp = sf_tprod+sp_tprod+sd_tprod+mf_tprod+mp_tprod+md_tprod+lp_tprod+ld_tprod;

%%
load([ipath 'Means_SSP585_empHP_2090-2100_' cfile '.mat'],...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod');

IStcb = sf_tmean+sp_tmean+sd_tmean+mf_tmean+mp_tmean+md_tmean+lp_tmean+ld_tmean;
IStcp = sf_tprod+sp_tprod+sd_tprod+mf_tprod+mp_tprod+md_tprod+lp_tprod+ld_tprod;

%% time
yH = 1950+(1/12):(1/12):2015;
yS = 2015+(1/12):(1/12):2101;

% ref period for Vianney is 1995-2014

tid = find(yH>1995 & yH<=2014);

GHbm = mean(GHtcb(tid));
GHpm = mean(GHtcp(tid));
IHbm = mean(IHtcb(tid));
IHpm = mean(IHtcp(tid));

Gtcp = [GHtcp GStcp];
Gtcb = [GHtcb GStcb];
Itcp = [IHtcp IStcp];
Itcb = [IHtcb IStcb];

pdGHtcb = (Gtcb-GHbm)./GHbm;
pdGHtcp = (Gtcp-GHpm)./GHpm;
pdIHtcb = (Itcb-IHbm)./IHbm;
pdIHtcp = (Itcp-IHpm)./IHpm;

yr = [yH yS];

%% percent diff from 1995-2014
figure(1)
subplot(2,2,1)
plot(yr,100*pdGHtcp,'r','LineWidth',1); hold on;
plot(yr,100*pdGHtcb,'b','LineWidth',1); hold on;
title('All fish consumers GFDL-FEISTY')
ylabel('% difference from 1995-2014')
legend('production','biomass')
xlim([1950 2100])

subplot(2,2,2)
plot(yr,100*pdIHtcp,'r','LineWidth',1); hold on;
plot(yr,100*pdIHtcb,'b','LineWidth',1); hold on;
title('All fish consumers IPSL-FEISTY')
ylabel('% difference from 1995-2014')
legend('production','biomass')
xlim([1950 2100])

stamp('')
print('-dpng',[ppath 'Hist_SSP585_ts_empHP_tcb_prod_pdiff_1994_2014.png'])


