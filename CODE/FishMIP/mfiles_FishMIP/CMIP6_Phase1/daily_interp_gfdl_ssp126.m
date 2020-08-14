% Make mat files of interpolated time series from GFDL
% SSP 126 2015-2100

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/ssp126/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-3
%tp: degC
%tb: degC

load([fpath 'gfdl_ssp126_temp100_monthly_2015_2100.mat'],'temp_100');
load([fpath 'gfdl_ssp126_temp_btm_monthly_2015_2100.mat'],'temp_btm');
load([fpath 'gfdl_ssp126_zmeso100_monthly_2015_2100.mat'],'zmeso_100');
load([fpath 'gfdl_ssp126_det_btm_monthly_2015_2100.mat']); %,'det_btm'

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
Tdays=1:365;
Time=Tdays(15:30:end);

%%
% index of water cells
[ni,nj,nt] = size(temp_btm);
WID = find(~isnan(temp_btm(:,:,1)));  % spatial index of water cells
NID = length(WID);                    % number of water cells 44564

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
        
        % meso zoo: from molC m-3 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % mult by 10 m depth interval for m-3 to m-2
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zm(j,:) = yi * 12.01 * 9.0 * 10;
        
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
    save([fpath 'Data_gfdl_ssp126_daily_',num2str(yr),'.mat'], 'ESM');
    
    
end




