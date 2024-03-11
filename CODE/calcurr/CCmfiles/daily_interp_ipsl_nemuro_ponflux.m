% Make mat files of interpolated time series from NEMURO CCE model
% using IPSL downscaled projections

clear 
close all

%%
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
%fpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/ESM_data/NEMURO/';
%spath = '/Users/cpetrik/Documents/NEMURO/';

load([fpath 'feisty_ipsl_gridspec.mat'],'TMO','TMO_bnds','tmo_units','BATHY')
load([fpath 'Data_grid_nemuro_ipsl.mat'],'GRD');

%%
time = (TMO/365)+1900;
yrs= 1980:2100;
Tdays=1:365;

mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

%% Units
%poc_bot: mmolN/m3/d %Might be mmolN/m2 now?
%lzoo: mmolN/m2
%pzoo: mmolN/m2
%tp: degC
%tb: degC


load([fpath 'feisty_ipsl_lzoo_1980-2100.mat'],'LZOO_INT_200M',...
    'lz_long_name','lz_units');

load([fpath 'feisty_ipsl_pon_1980-2100.mat'],'PON_BOT',...
    'pon_long_name','pon_units');

load([fpath 'feisty_ipsl_pzoo_1980-2100.mat'],'PZOO_INT_200M',...
    'pz_long_name','pz_units');

load([fpath 'feisty_ipsl_tb_1980-2100.mat'],'TEMP_BOT','tb_long_name','tb_units');

load([fpath 'feisty_ipsl_tp_1980-2100.mat'],'TEMP_AVG_200M',...
    'tp_long_name','tp_units','LAT','LON');

[ni,nj,nt] = size(TEMP_BOT);

%%
PON_BOT = -1.0 * PON_BOT; %was neg b/c downward flux

LZOO_INT_200M(LZOO_INT_200M<0) = 0.0;
PZOO_INT_200M(PZOO_INT_200M<0) = 0.0;
PON_BOT(PON_BOT<0) = 0.0;

%% Calc Martin curve for POC to bottom
% b = -0.863 for whole CCE in Martin et al. 1987

depth = repmat(BATHY,1,1,nt);
det_btm = PON_BOT .* ((depth/100).^-0.863);

%%
figure
pcolor(PON_BOT(:,:,1)); shading flat;
caxis([0 0.002])

figure
pcolor(det_btm(:,:,1)); shading flat;
caxis([0 0.002])
  
%%
WID = GRD.ID;  %1560 of these NaNs - maybe FW bodies on land?
NID = GRD.N;

for y = 1:length(yrs)
    yr = yrs(y)

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
    
    Tp = TEMP_AVG_200M(:,:,range);
    Tb = TEMP_BOT(:,:,range);
    Zm = LZOO_INT_200M(:,:,range);
    Zl = PZOO_INT_200M(:,:,range);
    det= det_btm(:,:,range);
    
    
    % setup FEISTY data files
    D_Tp  = zeros(NID,365);
    D_Tb  = zeros(NID,365);
    D_Zm  = zeros(NID,365);
    D_Zl  = zeros(NID,365);
    D_det = zeros(NID,365);
    
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
        
        % medium zoo: from mmolN m-2 to g(WW) m-2
        % 1e-3 mol in 1 mmol
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm(j,:) = yi * 1e-3 * (106.0/16.0) * 12.01 * 9.0;
        
        % large zoo: from mmolN m-2 to g(WW) m-2
        Y = squeeze(Zl(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zl(j,:) = yi * 1e-3 * (106.0/16.0) * 12.01 * 9.0;
        
        % detrital flux to benthos: from mmolN m-3 d-1 to g(WW) m-2 d-1
        % 1e-3 mol in 1 mmol
        % 106/16 mol C in 1 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % mult by 20 m depth interval for m-3 to m-2
        Y = squeeze(det(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_det(j,:) = yi * 1e-3  * (106.0/16.0) * 12.01 * 9.0 * 20;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_Zl(D_Zl<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
    ESM.Tp = D_Tp;
    ESM.Tb = D_Tb;
    ESM.Zm = D_Zm;
    ESM.Zl = D_Zl;
    ESM.det = D_det;
    
    % save
    save([spath 'Data_nemuro_ipsl_' num2str(yr) '.mat'], 'ESM');
    
end



