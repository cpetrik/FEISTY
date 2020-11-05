% check mesoz inputs

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

%% Pre

load([ppath 'ipsl_pi_temp_100_monthly_1950_2100.mat']); 

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;
temp_100 = fliplr(temp_100);
pTemp = temp_100;

pyr = yr(runs);

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
pT = temp_100(WID,:);

%% Compare to spinup
clear temp_100 runs yr

load([ppath 'ipsl_pi_temp_100_monthly_1850_1949.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;
temp_100 = fliplr(temp_100);
sTemp = temp_100;

syr = yr(spin);

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
sT = temp_100(WID,:);

%% Hist 100
clear temp_100 spin yr

load([hpath 'ipsl_hist_temp_100_monthly_1950_2014.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;
temp_100 = fliplr(temp_100);
hTemp = temp_100;

hyr = yr(runs);

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
hT = temp_100(WID,:);

%% SSP 126
clear temp_100 runs yr

load([fpath 'ipsl_ssp126_temp_100_monthly_2015_2100.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;
fTemp = temp_100;

fyr = yr;

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
fT = temp_100(WID,:);

%% SSP 585
clear temp_100 yr

load([spath 'ipsl_ssp585_temp_100_monthly_2015_2100.mat']);

temp_100 = double(temp_100);
temp_100(temp_100 > 1.0e19) = nan;
rTemp = temp_100;

ryr = yr;

[ni,nj,nt] = size(temp_100);
temp_100 = reshape(temp_100,ni*nj,nt);
rT = temp_100(WID,:);

%%
% load('ipsl_temp_ts.mat','syr','pyr','hyr','fyr','sT','pT','hT','fT','rT');

%% plot
figure(1)
plot(syr,nanmean(sT),'color',[0.5 0.5 0.5]); hold on
plot(pyr,nanmean(pT),'k'); hold on
plot(hyr,nanmean(hT),'b'); hold on
plot(fyr,nanmean(fT),'m'); hold on
plot(ryr,nanmean(rT),'r'); hold on
title('Temp')

%%
save('ipsl_temp_ts.mat','syr','pyr','hyr','fyr','sT','pT','hT','fT','rT');

%%
figure(2)
subplot(2,2,1)
pcolor(nanmean(pTemp,3)')
shading flat
cmocean('thermal')
colorbar
caxis([0 30])

subplot(2,2,2)
pcolor(nanmean(hTemp,3)')
shading flat
cmocean('thermal')
colorbar
caxis([0 30])

subplot(2,2,3)
pcolor(nanmean(fTemp,3)')
shading flat
cmocean('thermal')
colorbar
caxis([0 30])

subplot(2,2,4)
pcolor(nanmean(rTemp,3)')
shading flat
cmocean('thermal')
colorbar
caxis([0 30])

%%






