% check pre-indust inputs, because fish biomass too low

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/preindust/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC
% D_Zm = yi * 12.01 * 9.0;
% D_det = yi * 12.01 * 9.0 * 60 * 60 * 24;

load([fpath 'gfdl_pi_temp_100_monthly_1950_2100.mat'],'temp_100');
load([fpath 'gfdl_pi_temp_btm_monthly_1950_2100.mat'],'temp_btm');
load([fpath 'gfdl_pi_zmeso_100_monthly_1950_2100.mat'],'zmeso_100');
load([fpath 'gfdl_pi_det_btm_monthly_1950_2100.mat']); %,'det_btm'

temp_100 = double(temp_100);
temp_btm = double(temp_btm);
zmeso_100 = double(zmeso_100);
det_btm = double(det_btm);

%%
% temp_100(temp_100 > 1.0e19) = nan;
% temp_btm(temp_btm > 1.0e19) = nan;
% zmeso_100(zmeso_100 > 1.0e19) = nan;
% det_btm(det_btm > 1.0e19) = nan;

%%
ryr = yr(runs);

% index of water cells
[ni,nj,nt] = size(temp_100);
WID = find(~isnan(temp_100(:,:,1)));  % spatial index of water cells
NID = length(WID);

%%
temp_100 = reshape(temp_100,ni*nj,nt);
temp_btm = reshape(temp_btm,ni*nj,nt);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
det_btm = reshape(det_btm,ni*nj,nt);

pTp = temp_100(WID,:);
pTb = temp_btm(WID,:);
pZ = zmeso_100(WID,:);
pD = det_btm(WID,:);

%% Compare to spinup
clear temp_100 temp_btm zmeso_100 det_btm

load([fpath 'gfdl_pi_temp_100_monthly_1850_1949.mat'],'temp_100');
load([fpath 'gfdl_pi_temp_btm_monthly_1850_1949.mat'],'temp_btm');
load([fpath 'gfdl_pi_zmeso_100_monthly_1850_1949.mat'],'zmeso_100');
load([fpath 'gfdl_pi_det_btm_monthly_1850_1949.mat']); %,'det_btm'

temp_100 = double(temp_100);
temp_btm = double(temp_btm);
zmeso_100 = double(zmeso_100);
det_btm = double(det_btm);

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

%%
syr = yr(spin);

% index of water cells
[ni,nj,nt] = size(temp_100);
WID = find(~isnan(temp_100(:,:,1)));  % spatial index of water cells
NID = length(WID);

%%
temp_100 = reshape(temp_100,ni*nj,nt);
temp_btm = reshape(temp_btm,ni*nj,nt);
zmeso_100 = reshape(zmeso_100,ni*nj,nt);
det_btm = reshape(det_btm,ni*nj,nt);

sTp = temp_100(WID,:);
sTb = temp_btm(WID,:);
sZ = zmeso_100(WID,:);
sD = det_btm(WID,:);

%%
clear temp_100 temp_btm zmeso_100 det_btm

%% plot
figure(1)
plot(syr,nanmean(sTp),'k'); hold on
plot(ryr,nanmean(pTp),'b'); hold on
title('Ptemp')

figure(2)
plot(syr,nanmean(sTb),'k'); hold on
plot(ryr,nanmean(pTb),'b'); hold on
title('Btemp')

figure(3)
plot(syr,nanmean(sZ),'k'); hold on
plot(ryr,nanmean(pZ),'b'); hold on
title('Zoo')

figure(4)
plot(syr,nanmean(sD),'k'); hold on
plot(ryr,nanmean(pD),'b'); hold on
title('Det')




