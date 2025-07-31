% Make mat files of interpolated time series from GFDL
% Reanalysis-forced runs 1993-2019
% MOM6-NWA12

clear
close all

cpath='/project/Feisty/GCM_Data/MOM6-NWA12/';

%% Units
%poc flux: mol N m-2 s-1
%zoo: mol N m-2
%tp: degC
%tb: degC

load([cpath 'temp_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],'temp_100');
load([cpath 'btm_temp.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],'btm_temp');
load([cpath 'fntot_btm.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],'fntot_btm','fntot_btm_long_name','fntot_btm_units')
load([cpath 'nmdz_nlgz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat']);
load([cpath 'jhploss_mdz_lgz_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat']);

%%
temp_100(temp_100 > 1.0e19) = nan;
btm_temp(btm_temp > 1.0e19) = nan;
fntot_btm(fntot_btm > 1.0e19) = nan;
nmdz_100(nmdz_100 > 1.0e19) = nan;
nlgz_100(nlgz_100 > 1.0e19) = nan;
jhploss_nmdz_100(jhploss_nmdz_100 > 1.0e19) = nan;
jhploss_nlgz_100(jhploss_nlgz_100 > 1.0e19) = nan;

%%
load([fpath 'nwa_raw_ocean_static_gridspec.mat'],'geolon','geolat');
load([fpath 'Data_grid_mom6_nwa12.mat'],'GRD');

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

Tdays=1:365;

%%
% size
[ni,nj,nt] = size(fntot_btm);

% index of water cells
WID = GRD.ID;       % spatial index of water cells
NID = length(WID);  % number of water cells

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
    tb = (btm_temp(:,:,range));
    poc= (fntot_btm(:,:,range));
    zm = (nmdz_100(:,:,range));
    zl = (nlgz_100(:,:,range));
    hpm = (jhploss_nmdz_100(:,:,range));
    hpl = (jhploss_nlgz_100(:,:,range));

    % setup FEISTY data files
    Tp   = nan*zeros(NID,365);
    Tb   = nan*zeros(NID,365);
    det  = nan*zeros(NID,365);
    Zm   = nan*zeros(NID,365);
    Zl   = nan*zeros(NID,365);
    dZm  = nan*zeros(NID,365);
    dZl  = nan*zeros(NID,365);

    % interpolate to daily resolution
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

        % MZ: from molN m-2 to g(WW) m-2
        % 106 mol C to 16 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zm(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zm(j,:) = yi * (106/16) * 12.01  * 9.0;

        % LZ: from molN m-2 to g(WW) m-2
        % 106 mol C to 16 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zl(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zl(j,:) = yi * (106/16) * 12.01  * 9.0;

        % HPloss MZ: from molN m-2 s-1 to g(WW) m-2 d-1
        % 106 mol C to 16 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(hpm(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        dZm(j,:) = yi * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

        % HPloss LZ: from molN m-2 s-1 to g(WW) m-2 d-1
        % 106 mol C to 16 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(hpl(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        dZl(j,:) = yi * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

        % detrital flux to benthos: from molN m-2 s-1 to g(WW) m-2 d-1
        % 106 mol C to 16 mol N
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(poc(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        det(j,:) = yi * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

    end

    % Negative biomass or mortality loss from interp
    Zm(Zm<0)   = 0.0;
    Zl(Zl<0)   = 0.0;
    dZm(dZm<0) = 0.0;
    dZl(dZl<0) = 0.0;
    det(det<0) = 0.0;

    ESM.Tp  = Tp;
    ESM.Tb  = Tb;
    ESM.Zm  = Zm;
    ESM.Zl  = Zl;
    ESM.dZm = dZm;
    ESM.dZl = dZl;
    ESM.det = det;

    % save
    save([cpath 'Data_mom6_nwa12_daily_',num2str(YR),'.mat'], 'ESM','-v7.3');
    
end

%%

