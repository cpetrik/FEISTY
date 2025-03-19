% Make mat files of interpolated time series from CESM
% Hist 1950-2014
% 200m vertical integrations

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
ppath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/preindust/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'cesm2_hist_temp_btm_monthly_1850_2014.mat'],'temp_btm');
load([fpath 'cesm2_hist_temp_150_monthly_1850_2014.mat'],'temp_150');
load([fpath 'cesm2_hist_zooc_150_monthly_1850_2014.mat'],'zooc_150','units_vint');
load([fpath 'cesm2_hist_phyc_150_monthly_1850_2014.mat'],'phyc_150');
load([fpath 'cesm2_hist_diat_150_monthly_1850_2014.mat'],'diat_150');
load([fpath 'cesm2_hist_det_monthly_1850_2014.mat']); %,'det'

load([ppath 'cesm2-waccm_r1i1p1f1_picontrol_deptho_60arcmin_global_fx.mat'],'deptho')

%%
temp_btm(temp_btm > 1.0e19) = nan;
temp_150(temp_150 > 1.0e19) = nan;
zooc_150(zooc_150 > 1.0e19) = nan;
phyc_150(phyc_150 > 1.0e19) = nan;
diat_150(diat_150 > 1.0e19) = nan;
det(det > 1.0e19) = nan;

%% Calc zmeso from diat frac
Lfrac = double(diat_150) ./ double(phyc_150);
zmeso_150 = Lfrac .* double(zooc_150);

%% test that all same orientation
test1 = squeeze(double(temp_150(:,:,70)));
test2 = squeeze(double(temp_btm(:,:,70)));
test3 = squeeze(double(zmeso_150(:,:,70)));
test4 = squeeze(double(det(:,:,70)));

figure
subplot(2,2,1)
pcolor(test1); shading flat
subplot(2,2,2)
pcolor(test2); shading flat
subplot(2,2,3)
pcolor(test3); shading flat
subplot(2,2,4)
pcolor(test4); shading flat

figure
pcolor(deptho); shading flat

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
[ni,nj,nt] = size(det);
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

%%
Zoo = double(reshape(zooc_150,ni*nj,nt));
Phy = double(reshape(phyc_150,ni*nj,nt));
Dia = double(reshape(diat_150,ni*nj,nt));
frac = (reshape(Lfrac,ni*nj,nt));

Zoo = Zoo(WID,:);
Phy = Phy(WID,:);
Dia = Dia(WID,:);
frac = frac(WID,:);

hist_Zoo = mean(Zoo);
hist_Phy = mean(Phy);
hist_Dia = mean(Dia);
hist_frac = mean(frac);

figure
subplot(2,2,1)
plot(yr,hist_Phy,'r')

subplot(2,2,2)
plot(yr,hist_Dia,'b')

subplot(2,2,3)
plot(yr,hist_Zoo,'color',[0.75 0 0.5])

subplot(2,2,4)
plot(yr,hist_frac,'color',[0 0.5 0.75]) %goes up to 1e6, but should be constrained to 0-1

%% Lfrac
mfrac = mean(double(Lfrac),3);

figure
pcolor(mfrac); shading flat
%problem at edge of Antarctica
