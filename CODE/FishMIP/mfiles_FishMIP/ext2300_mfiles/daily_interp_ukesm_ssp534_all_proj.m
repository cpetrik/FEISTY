% Make mat files of interpolated time series from UKESM1-0-LL
% SSP 534-over 2101-2300
% 200 m vertical integrations
% bias-corrected btm temp

clear 
close all

fpath='/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'ukesm_ssp534_temp_btm_corrected_monthly_2040_2300.mat'],'temp_btm');

%%
load([fpath 'ukesm_ssp534_temp_200_monthly_2040_2100.mat'],'temp_200');
load([fpath 'ukesm_ssp534-over_zmeso_200_monthly_2040_2100.mat'],'zmeso_200');
load([fpath 'ukesm_ssp534_det_monthly_2040_2100.mat']);

tp1 = temp_200;
zm1 = zmeso_200;
dt1 = det;
time1 = time;
yr1 = yr;

clear zmeso_200 temp_200 det time yr

%%
load([fpath 'ukesm_ssp534_temp_200_monthly_2101_2300.mat'],'temp_200');
load([fpath 'ukesm_ssp534-over_zmeso_200_monthly_2101_2300.mat']);
load([fpath 'ukesm_ssp534_det_monthly_2101_2300.mat']);

tp2 = temp_200;
zm2 = zmeso_200;
dt2 = det;
time2 = time;
yr2 = yr;

clear zmeso_200 temp_200 det time yr

%%
[ni,nj,nt1] = size(tp1);
[~,~,nt] = size(temp_btm);

temp_200 = tp1;
temp_200(:,:,(nt1+1):nt)= tp2;

zmeso_200 = zm1;
zmeso_200(:,:,(nt1+1):nt)= zm2;

det = dt1;
det(:,:,(nt1+1):nt)= dt2;

time = [time1; time2];
yr = [yr1; yr2];

%%
temp_200(temp_200 > 1.0e19) = nan;
temp_btm(temp_btm > 1.0e19) = nan;
zmeso_200(zmeso_200 > 1.0e19) = nan;
det(det > 1.0e19) = nan;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):floor(yr(end));
Tdays=1:365;

% yrs = 1850:2014;
% Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(temp_200(:,:,70)));
test2 = squeeze(double(temp_btm(:,:,70)));
test3 = squeeze(double(zmeso_200(:,:,70)));
test4 = squeeze(double(det(:,:,70)));

% figure
% subplot(2,2,1)
% pcolor(test1); shading flat
% subplot(2,2,2)
% pcolor(test2); shading flat
% subplot(2,2,3)
% pcolor(test3); shading flat
% subplot(2,2,4)
% pcolor(test4); shading flat


%% index of water cells
%make GRD in another file later
% load('/project/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/Data_grid_ukesm.mat','GRD');
% WID = GRD.ID;
% NID = GRD.N;

% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));
NID = length(WID);

[ni,nj] = size(test4);

%%
for y = 1:nyrs
    ytime = yrs(y);

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
    
    Tp = double(temp_200(:,:,range));
    Tb = double(temp_btm(:,:,range));
    Zm = double(zmeso_200(:,:,range));
    Det= double(det(:,:,range));
    
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
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % meso zoo: from molC m-2 to g(WW) m-2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm(j,:) = yi * 12.01 * 9.0;
        
        % detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(Det(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
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
    save([fpath 'Data_ukesm_ssp534-over_daily_',num2str(ytime),'.mat'], 'ESM');
    
end

