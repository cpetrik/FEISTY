% check mesoz inputs

clear all
close all

ppath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/preindust/';
hpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/hist/';
fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp126/';
spath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/ssp585/';

%% Grid
load('/Volumes/MIP/Fish-MIP/CMIP6/GFDL/Data_grid_gfdl.mat','GRD');
NX = length(GRD.Z);
WID = GRD.ID;

%% Pre

load([ppath 'gfdl_pi_temp_100_monthly_1950_2100.mat']); 

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;

pyr = yr(runs);

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
pT = temp_100(WID,:);

%% Compare to spinup
clear temp_100 runs yr

load([ppath 'gfdl_pi_temp_100_monthly_1850_1949.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;

syr = yr(spin);

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
sT = temp_100(WID,:);

%% Hist 100
clear temp_100 spin yr

load([hpath 'gfdl_hist_temp_100_monthly_1950_2014.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;

hyr = yr(runs);

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
hT = temp_100(WID,:);

%% SSP 126
clear temp_100 runs yr

load([fpath 'gfdl_ssp126_temp_100_monthly_2015_2100.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;

fyr = yr;

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
fT = temp_100(WID,:);

%% SSP 585
clear temp_100 yr

load([spath 'gfdl_ssp585_temp_100_monthly_2015_2100.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;

ryr = yr;

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
rT = temp_100(WID,:);

%%
load('gfdl_temp_ts.mat','syr','pyr','hyr','fyr','sT','pT','hT','fT','rT');

%% plot
figure(1)
plot(syr,nanmean(sT),'color',[0.5 0.5 0.5]); hold on
plot(pyr,nanmean(pT),'k'); hold on
plot(hyr,nanmean(hT),'b'); hold on
plot(fyr,nanmean(fT),'m'); hold on
plot(ryr,nanmean(rT),'r'); hold on
title('Temp')

%%
save('gfdl_temp_ts.mat','syr','pyr','hyr','fyr','sT','pT','hT','fT','rT');

%%






