% Make mat files of interpolated time series from GFDL
% Preindust runs 1950-2100
% New vertical integrations

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');
load([fpath 'gridspec_POP_gx1v6.mat'],'mask');

%% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

%% Units
%poc flux: mmol/m^3 cm/s
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC
%tb: degC

% From nmolC cm-2 s-1 to g(WW) m-2 d-1
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

% From nmolC cm-2 s-1 to g(WW) m-2
% * 1e-9 * 1e4 * 12.01 * 9.0;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;
yrs = 1:nyrs;

Tdays=1:365;
Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(TEMP_150m(:,:,200)));
test2 = squeeze(double(TEMP_bottom(:,:,200)));
test3 = squeeze(double(LzooC_150m(:,:,200)));
test4 = squeeze(double(Lzoo_loss_150m(:,:,200)));
test5 = squeeze(double(POC_FLUX_IN_bottom(:,:,200)));

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

figure
subplot(2,3,1)
pcolor(test1)
shading flat
subplot(2,3,2)
pcolor(test2)
shading flat
subplot(2,3,3)
pcolor(test5)
shading flat
subplot(2,3,4)
pcolor(test3)
shading flat
subplot(2,3,5)
pcolor(test4)
shading flat
print('-dpng',[pp 'cesm_fosi.png'])

%%
for y = 1:nyrs
    yr = yrs(y)
    
    Tp  = (TEMP_150m(:,:,mstart(y):mend(y)));
    Tb  = (TEMP_bottom(:,:,mstart(y):mend(y)));
    Zm  = (LzooC_150m(:,:,mstart(y):mend(y)));
    dZm = (Lzoo_loss_150m(:,:,mstart(y):mend(y)));
    det = (POC_FLUX_IN_bottom(:,:,mstart(y):mend(y)));
    
    % index of water cells
    [ni,nj,nt] = size(TEMP_bottom);
    WID = find(~isnan(mask));  % spatial index of water cells
    NID = length(WID);         % number of water cells
    
    % setup FEISTY data files
    D_Tp  = nan*zeros(NID,365);
    D_Tb  = nan*zeros(NID,365);
    D_Zm  = nan*zeros(NID,365);
    D_dZm  = nan*zeros(NID,365);
    D_det = nan*zeros(NID,365);
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        % pelagic temperature (in Celcius)
        Y = squeeze(Tp(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % meso zoo: nmolC cm-2 to g(WW) m-2
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0;
        
        % medium zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(dZm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_dZm(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
        
        % detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(det(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_det(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_dZm(D_dZm<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    D_Zm(isnan(D_Zm)) = 0.0;
    D_dZm(isnan(D_dZm)) = 0.0;
    D_det(isnan(D_det)) = 0.0;
    
    ESM.Tp = D_Tp;
    ESM.Tb = D_Tb;
    ESM.Zm = D_Zm;
    ESM.dZm = D_dZm;
    ESM.det = D_det;
    
    % save
    save([fpath 'Data_cesm_fosi_daily_',num2str(yr),'.mat'], 'ESM');
    
    
end


