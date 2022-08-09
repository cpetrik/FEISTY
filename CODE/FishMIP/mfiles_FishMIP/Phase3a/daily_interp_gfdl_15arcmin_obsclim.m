% Make mat files of interpolated time series from GFDL
% Reanalysis-forced runs 1961-2010
% Obsclim

clear all
close all

fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'gfdl-mom6-cobalt2_obsclim_temp100_15arcmin_global_monthly_1961_2010.mat'],'temp_100');
load([fpath 'gfdl-mom6-cobalt2_obsclim_tob_15arcmin_global_monthly_1961_2010.mat'],'tob');
load([fpath 'gfdl-mom6-cobalt2_obsclim_zmeso100_15arcmin_global_monthly_1961_2010.mat'],'zmeso_100');
load([fpath 'gfdl-mom6-cobalt2_obsclim_expc-bot_15arcmin_global_monthly_1961_2010.mat']); %,'det_btm'

temp_100(temp_100 > 1.0e19) = nan;
tob(tob > 1.0e19) = nan;
zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

Tdays=1:365;

%% test that all same orientation
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';

figure
subplot(2,2,1)
pcolor(squeeze(double(temp_100(:,:,200)))); shading flat
subplot(2,2,2)
pcolor(squeeze(double(tob(:,:,200)))); shading flat
subplot(2,2,3)
pcolor(squeeze(double(zmeso_100(:,:,200)))); shading flat
subplot(2,2,4)
pcolor(squeeze(double(det_btm(:,:,200)))); shading flat
print('-dpng',[pp 'gfdl_obsclim_test.png'])

%%
% index of water cells
[ni,nj,nt] = size(det_btm);
WID = find(~isnan(det_btm(:,:,1)));  % spatial index of water cells
NID = length(WID);                    % number of water cells

%%
for y = 1%:nyrs
    YR = yrs(y)
    
    if y==1
        range = mstart(y):(mend(y)+1);
        Time=15:30:395;
    elseif y==nyrs
        range = (mstart(y)-1):mend(y);
        Time=-15:30:365;
    else
        range = (mstart(y)-1):(mend(y)+1);
        Time=-15:30:395;
    end
    
    Tp = (temp_100(:,:,range));
    Tb = (tob(:,:,range));
    Zm = (zmeso_100(:,:,range));
    det= (det_btm(:,:,range));
    
    % index of water cells
    [ni,nj,nt] = size(Tb);
    WID = find(~isnan(Tb(:,:,1)));  % spatial index of water cells
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
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % meso zoo: from molC m-2 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        D_Zm(j,:) = yi * 12.01 * 9.0;
        
        % detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(det(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
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
    save([fpath 'Data_gfdl_mom6_cobalt2_obsclim_15arcmin_daily_',num2str(YR),'.mat'], 'ESM','-v7.3');
    
    
end


