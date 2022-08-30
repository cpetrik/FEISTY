% Make mat files of interpolated time series from GFDL
% Reanalysis-forced runs 1961-2010
% Ctrlclim

clear all
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'gfdl-mom6-cobalt2_ctrlclim_temp100_15arcmin_global_monthly_1961_2010.mat'],'temp_100');
load([fpath 'gfdl-mom6-cobalt2_ctrlclim_tob_15arcmin_global_monthly_1961_2010.mat'],'tob');
load([fpath 'gfdl-mom6-cobalt2_ctrlclim_zmeso100_15arcmin_global_monthly_1961_2010.mat'],'zmeso_100');
load([fpath 'gfdl-mom6-cobalt2_ctrlclim_expc-bot_15arcmin_global_monthly_1961_2010.mat']); %,'det_btm'

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
print('-dpng',[pp 'gfdl_ctrlclim_test.png'])

%%
% index of water cells
[ni,nj,nt] = size(det_btm);
WID = find(~isnan(det_btm(:,:,1)));  % spatial index of water cells
NID = length(WID);                    % number of water cells

%%
for y = 1:nyrs
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
    
    tp = (temp_100(:,:,range));
    tb = (tob(:,:,range));
    zm = (zmeso_100(:,:,range));
    poc= (det_btm(:,:,range));
    
    % index of water cells
    [ni,nj,nt] = size(tb);
    WID = find(~isnan(tb(:,:,1)));  % spatial index of water cells
    NID = length(WID);              % number of water cells
    
    % setup FEISTY data files
    Tp  = nan*zeros(NID,365);
    Tb  = nan*zeros(NID,365);
    Zm  = nan*zeros(NID,365);
    det = nan*zeros(NID,365);
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        % pelagic temperature (in Celcius)
        Y = squeeze(tp(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(tb(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Tb(j,:) = yi;
        
        % meso zoo: from molC m-2 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zm(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zm(j,:) = yi * 12.01 * 9.0;
        
        % detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(poc(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        det(j,:) = yi * 12.01 * 9.0 * 60 * 60 * 24;
        
    end
    
    % Negative biomass or mortality loss from interp
    Zm(Zm<0) = 0.0;
    det(det<0) = 0.0;
    
    ESM.Tp = Tp;
    ESM.Tb = Tb;
    ESM.Zm = Zm;
    ESM.det = det;
    
    % save
    save([fpath 'Data_gfdl_mom6_cobalt2_ctrlclim_15arcmin_daily_',num2str(YR),'.mat'], 'ESM','-v7.3');
%     save([fpath 'gfdl_mom6_cobalt2_ctrlclim_15arcmin_daily_',num2str(YR),'.mat'],...
%         'Tp','Tb','Zm','det','-v7.3');
    
    
end


