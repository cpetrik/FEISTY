% Make mat files of interpolated time series from UKESM1-0-LL
% SSP 534-over 2101-2300
% 200 m vertical integrations

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'ukesm_ssp534_temp_btm_monthly_2040_2100.mat'],'temp_btm');
load([fpath 'ukesm_ssp534_temp_200_monthly_2040_2100.mat'],'temp_200');
load([fpath 'ukesm_ssp534_det_monthly_2040_2100.mat'])
% load([fpath 'ukesm_ssp534_zmeso_200_monthly_2040_2100.mat']) %DOESN'T EXIST

temp_200(temp_200 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
det(det > 1.0e19) = nan;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):floor(yr(end));
Tdays=1:365;

% yrs = 1850:2014;
% Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(temp_200(:,:,70)));
test2 = squeeze(double(temp_btm(:,:,70)));
%test3 = squeeze(double(zmeso_200(:,:,70)));
test4 = squeeze(double(det(:,:,70)));

figure
subplot(2,2,1)
pcolor(test1); shading flat
subplot(2,2,2)
pcolor(test2); shading flat
%subplot(2,2,3)
%pcolor(test3); shading flat
subplot(2,2,4)
pcolor(test4); shading flat


%% index of water cells
%make GRD in another file later
% load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/Data_grid_ukesm.mat','GRD');
% WID = GRD.ID;
% NID = GRD.N;

% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));
NID = length(WID);

[ni,nj] = size(test4);


%% Means over all grid cells
nt = length(time);

Tp = double(reshape(temp_200,ni*nj,nt));
Tb = double(reshape(temp_btm,ni*nj,nt));
% Zm = double(reshape(zmeso_200,ni*nj,nt));
Det= double(reshape(det,ni*nj,nt));

Tp = Tp(WID,:);
Tb = Tb(WID,:);
% Zm = Zm(WID,:);
Det= Det(WID,:);

ssp534_Tp = mean(Tp);
ssp534_Tb = mean(Tb);
% ssp534_Zm = mean(Zm);
ssp534_Det = mean(Det);

%%
figure
subplot(2,2,1)
plot(yr,ssp534_Tp,'r')

subplot(2,2,2)
plot(yr,ssp534_Tb,'b')

% subplot(2,2,3)
% plot(yr,ssp534_Zm,'color',[0.75 0 0.5])

subplot(2,2,4)
plot(yr,ssp534_Det,'color',[0 0.5 0.75])

%% save means
ssp534_yr = yr;
save([fpath 'Means_ukesm_ssp534_monthly_2040_2100.mat'], 'ssp534_Tp','ssp534_Tb',...
    'ssp534_Det','ssp534_yr'); %'ssp534_Zm',

