% Make mat files of interpolated time series from GFDL
% Reanalysis-forced runs 1950-2100
% Obsclim

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'gfdl_pi_temp_100_monthly_1950_2100.mat'],'temp_100');
load([fpath 'gfdl_pi_temp_btm_monthly_1950_2100.mat'],'temp_btm');
load([fpath 'gfdl_pi_zmeso_100_monthly_1950_2100.mat'],'zmeso_100');
load([fpath 'gfdl_pi_det_btm_monthly_1950_2100.mat']); %,'det_btm'

temp_100(temp_100 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

%%
mos = length(runs);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

ryr = yr(runs);
%yrs = floor(ryr(1)):ryr(end);
yrs = 1950:2100;

Tdays=1:365;
Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(temp_100(:,:,200)));
test2 = squeeze(double(temp_btm(:,:,200)));
test3 = squeeze(double(zmeso_100(:,:,200)));
test4 = squeeze(double(det_btm(:,:,200)));

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP6/';

figure
subplot(2,2,1)
pcolor(test1)
subplot(2,2,2)
pcolor(test2)
subplot(2,2,3)
pcolor(test3)
subplot(2,2,4)
pcolor(test4)
print('-dpng',[pp 'gfdl_pi_runs.png'])

%%
% index of water cells
[ni,nj,nt] = size(temp_100);
WID = find(~isnan(temp_100(:,:,1)));  % spatial index of water cells
NID = length(WID);                    % number of water cells

% setup FEISTY data files
% ESM.Tp  = nan*zeros(NID,365,nyrs);
% ESM.Tb  = nan*zeros(NID,365,nyrs);
% ESM.Zm  = nan*zeros(NID,365,nyrs);
% ESM.det = nan*zeros(NID,365,nyrs);

%%
for y = 1:nyrs
    yr = yrs(y)
    
    Tp = double(temp_100(:,:,mstart(y):mend(y)));
    Tb = double(temp_btm(:,:,mstart(y):mend(y)));
    Zm = double(zmeso_100(:,:,mstart(y):mend(y)));
    det= double(det_btm(:,:,mstart(y):mend(y)));
    
    % index of water cells
    [ni,nj,nt] = size(Tp);
    WID = find(~isnan(Tp(:,:,1)));  % spatial index of water cells
    NID = length(WID);              % number of water cells
    
    % setup FEISTY data files
    D_Tp  = nan*zeros(NID,365);
    D_Tb  = nan*zeros(NID,365);
    D_Zm  = nan*zeros(NID,365);
    D_det = nan*zeros(NID,365);
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        % pelagic temperature (in Celcius)
        Y = squeeze(Tp(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % meso zoo: from molC m-2 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zm(j,:) = yi * 12.01 * 9.0;
        
        % detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(det(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_det(j,:) = yi * 12.01 * 9.0 * 60 * 60 * 24;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    ESM.Tp = D_Tp;
    ESM.Tb = D_Tb;
    ESM.Zm = D_Zm;
    ESM.det = D_det;
    
    % save
    save([fpath 'Data_gfdl_pi_daily_',num2str(yr),'.mat'], 'ESM');
    
    
end


