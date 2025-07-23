% Bias correct CESM-WACCM Historic zooc
% Zmeso from diatom frac of zooc
% Bias-corrected with 50-yr mean (1965-2014) against obsGLMM

clear
close all

%% Hist
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
load([wpath 'gridspec_cesm2_cmip6_2300.mat']);
load([wpath 'Data_grid_cesm2_cmip6_2300.mat']);

load([wpath 'hist/cesm2_hist_biascorr_zooc_zmeso_diatfraconly_molC_monthly_1850_2014.mat'],...
    'diffZmeso');

%% SSP 126
load([wpath 'ssp126/cesm2_ssp126_zooc_150_monthly_2015_2299.mat']);
load([wpath 'ssp126/cesm2_ssp126_phyc_150_monthly_2015_2299.mat'],'phyc_150');
load([wpath 'ssp126/cesm2_ssp126_diat_150_monthly_2015_2299.mat'],'diat_150');

zooc_150(zooc_150 > 1.0e19) = nan;
phyc_150(phyc_150 > 1.0e19) = nan;
diat_150(diat_150 > 1.0e19) = nan;

% Calc zmeso from diat frac
Lfrac = diat_150 ./ phyc_150;
Lfrac(Lfrac>1) = 1.0;
Lfrac(Lfrac<0) = 0.0;

zmeso_150 = (Lfrac .* double(zooc_150));

%% bias-correct
[~,~,nt] = size(zmeso_150);
zmeso_corr = zmeso_150 - repmat(diffZmeso,1,1,nt);

%save
save([wpath 'ssp126/cesm2_ssp126_zmeso150_biascorr_monthly_2015_2299.mat'],'zmeso_150',...
    'Lfrac','zmeso_corr','units_vint','yr')

clear zmeso_150 zmeso_corr

%% SSP585
load([wpath 'ssp585/cesm2_ssp585_zooc_150_monthly_2015_2299.mat']);
load([wpath 'ssp585/cesm2_ssp585_phyc_150_monthly_2015_2299.mat'],'phyc_150');
load([wpath 'ssp585/cesm2_ssp585_diat_150_monthly_2015_2299.mat'],'diat_150');

%
zooc_150(zooc_150 > 1.0e19) = nan;
phyc_150(phyc_150 > 1.0e19) = nan;
diat_150(diat_150 > 1.0e19) = nan;

% Calc zmeso from diat frac
Lfrac = diat_150 ./ phyc_150;
Lfrac(Lfrac>1) = 1.0;
Lfrac(Lfrac<0) = 0.0;

zmeso_150 = (Lfrac .* double(zooc_150));

%% bias-correct
[~,~,nt] = size(zmeso_150);
zmeso_corr = zmeso_150 - repmat(diffZmeso,1,1,nt);

%save
save([wpath 'ssp585/cesm2_ssp585_zmeso150_biascorr_monthly_2015_2299.mat'],'zmeso_150',...
    'Lfrac','zmeso_corr','units_vint','yr')

clear zmeso_150 zmeso_corr

%% SSP534
load([wpath 'ssp534over/cesm2_ssp534-over_phyc_150_monthly_2040_2299.mat'],'phyc_150');
load([wpath 'ssp534over/cesm2_ssp534-over_diat_150_monthly_2040_2299.mat'],'diat_150');
load([wpath 'ssp534over/cesm2_ssp534-over_zooc_150_monthly_2040_2299.mat']); %,'zooc_150','units_vint');

zooc_150(zooc_150 > 1.0e19) = nan;
phyc_150(phyc_150 > 1.0e19) = nan;
diat_150(diat_150 > 1.0e19) = nan;

%calc zmeso from diat frac
Lfrac = diat_150 ./ phyc_150;
Lfrac(Lfrac>1) = 1.0;
Lfrac(Lfrac<0) = 0.0;

zmeso_150 = (Lfrac .* double(zooc_150));

%% bias-correct
[~,~,nt] = size(zmeso_150);
zmeso_corr = zmeso_150 - repmat(diffZmeso,1,1,nt);

%save
save([wpath 'ssp534over/cesm2_ssp534-over_zmeso150_biascorr_monthly_2040_2299.mat'],'zmeso_150',...
    'Lfrac','zmeso_corr','units_vint','yr')

clear zmeso_150 zmeso_corr
