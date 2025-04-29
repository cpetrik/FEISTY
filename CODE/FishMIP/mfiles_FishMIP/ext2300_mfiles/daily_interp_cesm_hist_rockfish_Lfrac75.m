% Make mat files of interpolated time series from CESM
% Hist 1950-2014
% 200m vertical integrations

clear 
close all

fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
ppath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/preindust/';

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
Lfrac(Lfrac>1) = 1.0;
Lfrac(Lfrac<0) = 0.0;

zmeso_150 = (Lfrac .* double(zooc_150)) + (0.75 * (1-Lfrac) .* double(zooc_150));

save([fpath 'cesm2_hist_zmeso_150_Lfrac75_monthly_1850_2014.mat'],'zmeso_150',...
    'Lfrac','units_vint','time','yr')

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
test1 = squeeze(double(temp_150(:,:,70)));
test2 = squeeze(double(temp_btm(:,:,70)));
test3 = squeeze(double(zmeso_150(:,:,70)));
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
% 
% figure
% pcolor(deptho); shading flat

%% index of water cells
%make GRD in another file later
% load('/project/Feisty/Fish-MIP/CMIP6/IPSL/Data_grid_cesm.mat','GRD');
% WID = GRD.ID;
% NID = GRD.N;

% Use btm det to find ocean cells
WID = find(~isnan(test4(:)));
NID = length(WID);

[ni,nj] = size(test4);

%%
for y = 100 %:nyrs
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

    Tp = double(temp_150(:,:,range));
    Tb = double(temp_btm(:,:,range));
    Zm = double(zmeso_150(:,:,range));
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
    save([fpath 'Data_cesm_hist_daily_',num2str(ytime),'_Lfrac75.mat'], 'ESM');


end



