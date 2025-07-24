% Make mat files of interpolated time series from CESM2-WACCM
% SSP 534-over 2101-2299
% 200 m vertical integrations
% bias-corrected btm temp
% Zmeso from diatom frac of zooc
% Bias-corrected with 50-yr mean (1965-2014) against obsGLMM

clear 
close all

%fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';
fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'cesm2_ssp534-over_temp_btm_corrected_monthly_2040_2299.mat'],'temp_btm');
load([fpath 'cesm2_ssp534-over_temp_150_monthly_2040_2299.mat'],'temp_150');
load([fpath 'cesm2_ssp534-over_det_monthly_2040_2299.mat'],'det')
load([fpath 'cesm2_ssp534-over_zooc_monthly_2040_2299.mat']); %,'zooc_150','units_vint');

temp_150(temp_150 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zooc_150(zooc_150 > 1.0e19) = nan;
det(det > 1.0e19) = nan;

%%
mos = length(yr);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):floor(yr(end));
Tdays=1:365;

% yrs = 1850:2014;
% Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(temp_150(:,:,70)));
test2 = squeeze(double(temp_btm(:,:,70)));
test3 = squeeze(double(zooc_150(:,:,70)));
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

%% index of water cells
%make GRD in another file later
% load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/Data_grid_cesm.mat','GRD');
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
Zm = double(reshape(zooc_150,ni*nj,nt));
Det= double(reshape(det,ni*nj,nt));

Tp = Tp(WID,:);
Tb = Tb(WID,:);
Zm = Zm(WID,:);
Det= Det(WID,:);

ssp534_Tp = mean(Tp);
ssp534_Tb = mean(Tb);
ssp534_Zm = mean(Zm,'omitnan');
ssp534_Det = mean(Det);

%%
figure
subplot(2,2,1)
plot(yr,ssp534_Tp,'r')

subplot(2,2,2)
plot(yr,ssp534_Tb,'b')

subplot(2,2,3)
plot(yr,ssp534_Zm,'color',[0.75 0 0.5])

subplot(2,2,4)
plot(yr,ssp534_Det,'color',[0 0.5 0.75])

%% save means
ssp534_yr = yr;
save([fpath 'Means_cesm2_ssp534_zooc_monthly_2040_2299.mat'], 'ssp534_Tp','ssp534_Tb',...
    'ssp534_Zm','ssp534_Det','ssp534_yr');
