% Make mat files of interpolated time series from CORE-forced ESM2M
% Hist 1950-2007

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%% Units
%det btm: mol N m-2 s-1
%zoo: integrated to mol N m-2  
%zoo loss: mol N m-2 s-1
%tp: degC
%tb: degC

load([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100');
load([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb');
load([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100');
load([fpath 'ocean_cobalt_lz100_monthly_1950_2007.mat'],'lz_100');
load([fpath 'ocean_cobalt_hploss_mz100_monthly_1950_2007.mat'],'hploss_mz_100');
load([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100');
load([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat']); %,'det_btm'

%%
mos = length(runs);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

ryr = yr(runs);
%yrs = ceil(ryr(1)):round(ryr(end));
yrs = 1950:2007;
Tdays=1:365;
Time=Tdays(15:30:end);

%% index of water cells
load([fpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
load([fpath 'ocean_cobalt_grid.mat'],'ht')

WID = GRD.ID;
NID = GRD.N;
[ni,nj] = size(ht);

%%
for y = 1:nyrs
    yr = yrs(y)
    
    Tp = double(tp_100(:,:,mstart(y):mend(y)));
    Tb = double(tb(:,:,mstart(y):mend(y)));
    Zm = double(mz_100(:,:,mstart(y):mend(y)));
    Zl = double(lz_100(:,:,mstart(y):mend(y)));
    dZm = double(hploss_mz_100(:,:,mstart(y):mend(y)));
    dZl = double(hploss_lz_100(:,:,mstart(y):mend(y)));
    det= double(det_btm(:,:,mstart(y):mend(y)));
    
    % setup FEISTY data files
    D_Tp  = nan*zeros(NID,365);
    D_Tb  = nan*zeros(NID,365);
    D_Zm  = nan*zeros(NID,365);
    D_Zl  = nan*zeros(NID,365);
    D_dZm  = nan*zeros(NID,365);
    D_dZl  = nan*zeros(NID,365);
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
        
        % medium zoo: from mol N m-2 to g(WW) m-2
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zm(j,:) = yi * (106.0/16.0) * 12.01 * 9.0;
        
        % large zoo: from mol N m-2 to g(WW) m-2
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zl(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_Zl(j,:) = yi * (106.0/16.0) * 12.01 * 9.0;
        
        % medium zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(dZm(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_dZm(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
        
        % large zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(dZl(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_dZl(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
        
        % detrital flux to benthos: from mol C m-2 s-1 to g(WW) m-2 d-1
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(det(m,n,:));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_det(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_Zl(D_Zl<0) = 0.0;
    D_dZm(D_dZm<0) = 0.0;
    D_dZl(D_dZl<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    COBALT.Tp = D_Tp;
    COBALT.Tb = D_Tb;
    COBALT.Zm = D_Zm;
    COBALT.Zl = D_Zl;
    COBALT.dZm = D_dZm;
    COBALT.dZl = D_dZl;
    COBALT.det = D_det;
    
    % save
    save([fpath 'Data_ocean_cobalt_daily_',num2str(yr),'.mat'], 'COBALT');
    
end



