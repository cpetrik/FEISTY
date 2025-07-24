% Make mat files of interpolated time series from UKESM1-0-LL
% Hist 1950-2014
% 200m vertical integrations

clear 
close all

fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'ukesm_hist_temp_btm_monthly_1850_2014.mat'],'temp_btm');
load([fpath 'ukesm_hist_temp_200_monthly_1850_2014.mat'],'temp_200');
load([fpath 'ukesm_hist_zmeso_200_monthly_1850_2014.mat'],'zmeso_200','units_vint');
load([fpath 'ukesm_hist_det_monthly_1850_2014.mat']); %,'det'

load('/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ukesm_isimip_depth.mat','deptho')

%%
temp_200(temp_200 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_200(zmeso_200 > 1.0e19) = nan;
det(det > 1.0e19) = nan;
deptho(deptho > 1.0e19) = nan;

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
test3 = squeeze(double(zmeso_200(:,:,70)));
test4 = squeeze(double(det(:,:,70)));

% figure
% subplot(2,2,1)
% pcolor(test1); shading flat
% subplot(2,2,2)
% pcolor(test2); shading flat
% subplot(2,2,3)
% pcolor(test3); shading flat
% subplot(2,2,4)
% pcolor(test4); shading flat
% 
% figure
% pcolor(deptho); shading flat

%% flip depth
depth = fliplr(deptho);

% figure
% pcolor(depth); shading flat

%% index of water cells
%make GRD in another file later
% load('/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/Data_grid_ukesm.mat','GRD');
% WID = GRD.ID;
% NID = GRD.N;

% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));
NID = length(WID);

[ni,nj] = size(test4);


%% Means over all grid cells
nt = length(yr);

Tp = double(reshape(temp_150,ni*nj,nt));
Tb = double(reshape(temp_btm,ni*nj,nt));
Zm = double(reshape(zmeso_150,ni*nj,nt));
Det= double(reshape(det,ni*nj,nt));

Tp = Tp(WID,:);
Tb = Tb(WID,:);
Zm = Zm(WID,:);
Det= Det(WID,:);

hist_Tp = mean(Tp);
hist_Tb = mean(Tb);
hist_Zm = mean(Zm);
hist_Det = mean(Det);

%%
figure
subplot(2,2,1)
plot(yr,hist_Tp,'r')

subplot(2,2,2)
plot(yr,hist_Tb,'b')

subplot(2,2,3)
plot(yr,hist_Zm,'color',[0.75 0 0.5])

subplot(2,2,4)
plot(yr,hist_Det,'color',[0 0.5 0.75])

%% save means
hist_yr = yr;
save([fpath 'Means_ukesm_hist_monthly_1850_2014.mat'], 'hist_Tp','hist_Tb',...
    'hist_Zm','hist_Det','hist_yr');


