% Make mat files of interpolated time series from GFDL
% Historic CMIP6 scenario 1880-2014
% Use MZ and LZ instead of one mesozoo

clear
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/hist/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol N m-2 and mol N m-2 s-1
%tp: degC
%tb: degC

load([fpath 'gfdl_hist_temp_100_monthly_1950_2014.mat'],'temp_100');
load([fpath 'gfdl_hist_temp_btm_monthly_1950_2014.mat'],'temp_btm');
load([fpath 'gfdl_hist_det_btm_monthly_1950_2014.mat']);

load([fpath 'gfdl_hist_int100_2zmeso_reorient_monthly_1950_2014.mat']);

temp_100 = double(temp_100);
temp_btm = double(temp_btm);
det_btm = double(det_btm);

%%
gpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/GFDL/';
load([gpath 'Data_grid_gfdl.mat'],'GRD');
load([gpath 'gridspec_gfdl_cmip6.mat']);

%%
yr65 = find(year>=1950);
time = time(yr65);
yr = year(yr65);

mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

Tdays=1:365;

%% test that all same orientation
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/FishMIP/FishMIP6/';

figure
subplot(3,3,1)
pcolor(squeeze((temp_100(:,:,200)))); shading flat
clim([0 30])
title('TP')

subplot(3,3,2)
pcolor(squeeze((temp_btm(:,:,200)))); shading flat
clim([0 30])
title('TB')

subplot(3,3,3)
pcolor(squeeze((mdz_100(:,:,200)))* (106/16) * 12.01 * 9.0); shading flat
clim([0 10])
title('MZ')

subplot(3,3,4)
pcolor(squeeze((lgz_100(:,:,200)))* (106/16) * 12.01 * 9.0); shading flat
clim([0 10])
title('LZ')

subplot(3,3,5)
pcolor(squeeze((hp_mdz_100(:,:,200)))* (106/16) * 12.01 * 9.0 * 60 * 60 * 24); shading flat
clim([0 1])
title('hpM')

subplot(3,3,6)
pcolor(squeeze((hp_lgz_100(:,:,200)))* (106/16) * 12.01 * 9.0 * 60 * 60 * 24); shading flat
clim([0 1])
title('hpL')

subplot(3,3,7)
pcolor(squeeze((det_btm(:,:,200)))* 12.01 * 9.0 * 60 * 60 * 24); shading flat
clim([0 5])
title('Det')

print('-dpng',[pp 'gfdl_cmip6_hist_2meso_test.png'])

%%
% index of water cells
[ni,nj,nt] = size(det_btm);
WID = find(~isnan(det_btm(:,:,1)));  % spatial index of water cells
NID = length(WID);                   % number of water cells = 44564

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
    tb = (temp_btm(:,:,range));
    poc= (det_btm(:,:,range));
    zm = (mdz_100(:,:,range));
    zl = (lgz_100(:,:,range));
    hpm = (hp_mdz_100(:,:,range));
    hpl = (hp_lgz_100(:,:,range));

   
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
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zm(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zm(j,:) = yi * (106/16) * 12.01 * 9.0;

        % LZ: from molN m-2 to g(WW) m-2
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zl(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zl(j,:) = yi * (106/16) * 12.01 * 9.0;

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

        % detrital flux to benthos: from molC m-2 s-1 to g(WW) m-2 d-1
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        % 60*60*24 sec in a day
        Y = squeeze(poc(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        det(j,:) = yi * 12.01 * 9.0 * 60 * 60 * 24;

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
    save([fpath 'Data_gfdl_cmip6_hist_2meso_daily_',num2str(YR),'.mat'], 'ESM','-v7.3');
    

end


