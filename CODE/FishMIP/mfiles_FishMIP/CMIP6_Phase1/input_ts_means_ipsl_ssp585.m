% Make mat files of time series from IPSL inputs
% SSP 585 2015-2100
% New vertical integrations

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/ssp585/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'ipsl_ssp585_temp_100_monthly_2015_2100.mat'],'temp_100');
load([fpath 'ipsl_ssp585_temp_btm_monthly_2015_2100.mat'],'temp_btm');
load([fpath 'ipsl_ssp585_zmeso_100_monthly_2015_2100.mat'],'zmeso_100');
load([fpath 'ipsl_ssp585_det_btm_monthly_2015_2100.mat']); %,'det_btm'

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;

nyrs = mos/12;
yrs = 2015:2100;

%% flip
temp_btm = fliplr(temp_btm);
det_btm = fliplr(det_btm);

%% test that all same orientation
test1 = squeeze(double(temp_100(:,:,120)));
test2 = squeeze(double(temp_btm(:,:,120)));
test3 = squeeze(double(zmeso_100(:,:,120)));
test4 = squeeze(double(det_btm(:,:,120)));

figure
subplot(2,2,1)
pcolor(test1)
shading flat
subplot(2,2,2)
pcolor(test2)
shading flat
subplot(2,2,3)
pcolor(test3)
shading flat
subplot(2,2,4)
pcolor(test4)
shading flat

%%
Tp = double(temp_100);
Tb = double(temp_btm);
Zm = double(zmeso_100);
det= double(det_btm);

% index of water cells
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/Data_grid_ipsl.mat','GRD');
load('/Volumes/MIP/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat','deptho')

WID = GRD.ID;
NID = GRD.N;
[ni,nj] = size(deptho);
[ti,tj,nt] = size(Tp);

%% Convert units
% meso zoo: from molC m-2 to g(WW) m-2
D_Zm = Zm * 12.01 * 9.0;

% detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
D_det = det * 12.01 * 9.0 * 60 * 60 * 24;

% Negative biomass or mortality loss from interp
D_Zm(D_Zm<0) = 0.0;
D_det(D_det<0) = 0.0;

%% time series
Tp = reshape(Tp,ni*nj,nt);
Tb = reshape(Tb,ni*nj,nt);
det = reshape(D_det,ni*nj,nt);
mzp = reshape(D_Zm,ni*nj,nt);

% Every 1 years
st=1:12:nt;
en=12:12:nt;
for n=1:length(st)
    Tp_1yr_ssp585(:,n)=nanmean(Tp(WID,st(n):en(n)),2);
    Tb_1yr_ssp585(:,n)=nanmean(Tb(WID,st(n):en(n)),2);
    det_1yr_ssp585(:,n)=nanmean(det(WID,st(n):en(n)),2);
    mz_1yr_ssp585(:,n)=nanmean(mzp(WID,st(n):en(n)),2);
end

%% save
s5time=time;
s5yrs=yrs;
save([fpath 'ipsl_ssp585_inputs_ts_nanmeans.mat'],'s5time','s5yrs',...
    'Tp_1yr_ssp585','Tb_1yr_ssp585','det_1yr_ssp585','mz_1yr_ssp585');



