% check mesoz inputs, because hist fish biomass very high

clear all
close all

ppath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/preindust/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp126/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';

%% Grid
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
NX = length(GRD.Z);
WID = GRD.ID;

%% Units
%zoo: mol C m-2
% D_Zm = yi * 12.01 * 9.0;

load([ppath 'ipsl_pi_zmeso_100_monthly_1950_2100.mat']); 

zmeso_100 = double(zmeso_100);
zmeso_100(zmeso_100 > 1.0e19) = nan;

pyr = yr(runs);

[ni,nj,nt] = size(zmeso_100);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
pZ = zmeso_100(WID,:);

%% Compare to spinup
clear zmeso_100 runs yr

load([ppath 'ipsl_pi_zmeso_100_monthly_1850_1949.mat']);

zmeso_100 = double(zmeso_100);
zmeso_100(zmeso_100 > 1.0e19) = nan;

syr = yr(spin);

[ni,nj,nt] = size(zmeso_100);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
sZ = zmeso_100(WID,:);

%% Hist 100
clear zmeso_100 spin yr

load([hpath 'ipsl_hist_zmeso_100_monthly_1950_2014.mat']);

zmeso_100 = double(zmeso_100);
zmeso_100(zmeso_100 > 1.0e19) = nan;

hyr = yr(runs);

[ni,nj,nt] = size(zmeso_100);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
hZ = zmeso_100(WID,:);

%% Hist vint
clear zmeso_100 yr

load([hpath 'ipsl_hist_zmeso_vint_monthly_1950_2014.mat']);

zmeso_vint = double(zmeso_vint);
zmeso_vint(zmeso_vint > 1.0e19) = nan;

[ni,nj,nt] = size(zmeso_vint);
zmeso_vint = reshape(zmeso_vint,ni*nj,nt);
hvZ = zmeso_vint(WID,:);

%% SSP 126
clear zmeso_100 runs yr

load([fpath 'ipsl_ssp126_zmeso_100_monthly_2015_2100.mat']);

zmeso_100 = double(zmeso_100);
zmeso_100(zmeso_100 > 1.0e19) = nan;

fyr = yr;

[ni,nj,nt] = size(zmeso_100);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
fZ = zmeso_100(WID,:);

%% SSP 585
clear zmeso_100 yr

load([spath 'ipsl_ssp585_zmeso_100_monthly_2015_2100.mat']);

zmeso_100 = double(zmeso_100);
zmeso_100(zmeso_100 > 1.0e19) = nan;

ryr = yr;

[ni,nj,nt] = size(zmeso_100);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
rZ = zmeso_100(WID,:);

%% plot
figure(1)
plot(syr,nanmean(sZ),'color',[0.5 0.5 0.5]); hold on
plot(pyr,nanmean(pZ),'k'); hold on
plot(hyr,nanmean(hZ),'b'); hold on
%plot(hyr,nanmean(hvZ),'c'); hold on
plot(fyr,nanmean(fZ),'m'); hold on
plot(ryr,nanmean(rZ),'r'); hold on
title('Zoo')





