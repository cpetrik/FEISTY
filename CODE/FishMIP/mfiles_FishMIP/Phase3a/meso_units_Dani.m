% Test values given to Dani

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

load([fpath 'gfdl-mom6-cobalt2_15arcmin_mdz_lgz_100m_global_monthly_1961_2010.mat'])
load([fpath 'gfdl-mom6_cobalt2_15arcmin_HPloss_mdz_lgz_100m_month_1961_2010.mat'])

%load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
%load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])

%% take means of first year and compare to daily interp data
MZ = squeeze(mean(nmdz_100(:,:,1:12),1,'omitnan'));
LZ = squeeze(mean(nlgz_100(:,:,1:12),1,'omitnan'));
hpMZ = squeeze(mean(hploss_nmdz_100(:,:,1:12),1,'omitnan'));
hpLZ = squeeze(mean(hploss_nlgz_100(:,:,1:12),1,'omitnan'));

MZ = mean(MZ,1,'omitnan');
LZ = mean(LZ,1,'omitnan');
hpMZ = mean(hpMZ,1,'omitnan');
hpLZ = mean(hpLZ,1,'omitnan');

%% units
%from gC m-2 to g(WW) m-2
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
MZ = MZ * 9.0;
LZ = LZ * 9.0;

%from gC m-2 s-1 to g(WW) m-2 d-1
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% 60*60*24 sec in a day
hpMZ = hpMZ * 9.0 * 60 * 60 * 24 *1e-3; %./ (1/1e-3);
hpLZ = hpLZ * 9.0 * 60 * 60 * 24 *1e-3; %./ (1/1e-3);

%%
load([cpath 'gfdl-mom6-cobalt2_obsclim_zmeso100_15arcmin_global_monthly_1961_2010.mat']);

%%
load([cpath 'Data_gfdl_mom6_cobalt2_obsclim_15arcmin_daily_1961.mat']);

% HP loss is an empirical fitted fn of biomass and temp
MZloss = 10 .^ (-2.925 + 1.964.*log10(ESM.Zm+eps) + 1.958e-2.*ESM.Tp);

totZ = mean(ESM.Zm,1,'omitnan');
hpTZ = mean(MZloss,1,'omitnan');

%% plot
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
cmo = cumsum(MNTH);
cmo = cmo-15;

figure(1) %7-9
plot(cmo, MZ+LZ, 'b'); hold on
plot(1:365,totZ,'r')


figure(2)
plot(cmo, hpMZ+hpLZ, 'b'); hold on %100-190 = 100x too high
plot(1:365,hpTZ,'r') %0.14-0.22
