% Make mat files of interpolated time series from GFDL
% Reanalysis-forced runs 1961-2010
% Obsclim
% Use MZ and LZ instead of one mesozoo

clear
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';

%% Units
%poc flux: mmol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

load([fpath 'gfdl-mom6-cobalt2_obsclim_temp100_15arcmin_global_monthly_1961_2010.mat'],'temp_100');
load([fpath 'gfdl-mom6-cobalt2_obsclim_tob_15arcmin_global_monthly_1961_2010.mat'],'tob');
%load([fpath 'gfdl-mom6-cobalt2_obsclim_zmeso100_15arcmin_global_monthly_1961_2010.mat'],'zmeso_100');
load([fpath 'gfdl-mom6-cobalt2_obsclim_expc-bot_15arcmin_global_monthly_1961_2010.mat'],...
    'det_btm','yr','lat','lon','LAT','LON')

temp_100(temp_100 > 1.0e19) = nan;
tob(tob > 1.0e19) = nan;
%zmeso_100(zmeso_100 > 1.0e19) = nan;
det_btm(det_btm > 1.0e19) = nan;

%%
gLAT = LAT;
gLON = LON;

%% 2 meso
load([cpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.mat']);  %units = 'molN/m2/s'
load([cpath '19610101-20101231.ocean_cobalt_tracers_int100_FishMIP_remapped.mat']); %units = 'gC/m2'
%nan = 1.000000020040877e+20

nmdz_100(nmdz_100 > 1.0e19) = nan;
nlgz_100(nlgz_100 > 1.0e19) = nan;
hploss_nmdz_100(hploss_nmdz_100 > 1.0e19) = nan;
hploss_nlgz_100(hploss_nlgz_100 > 1.0e19) = nan;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

Tdays=1:365;

%% test that all same orientation
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/FishMIP/Phase3a/';

figure
subplot(3,3,1)
pcolor(squeeze(double(temp_100(:,:,200)))); shading flat
clim([0 30])
title('TP')

subplot(3,3,2)
pcolor(squeeze(double(tob(:,:,200)))); shading flat
clim([0 30])
title('TB')

subplot(3,3,3)
pcolor(squeeze(double(nmdz_100(:,:,200)))* 9.0); shading flat
clim([0 10])
title('MZ')

subplot(3,3,4)
pcolor(squeeze(double(nlgz_100(:,:,200)))* 9.0); shading flat
clim([0 10])
title('LZ')

subplot(3,3,5)
pcolor(squeeze(double(hploss_nmdz_100(:,:,200)))* (106/16) * 12.01 * 9.0 * 60 * 60 * 24); shading flat
clim([0 1])
title('hpM')

subplot(3,3,6)
pcolor(squeeze(double(hploss_nlgz_100(:,:,200)))* (106/16) * 12.01 * 9.0 * 60 * 60 * 24); shading flat
clim([0 1])
title('hpL')

subplot(3,3,7)
pcolor(squeeze(double(det_btm(:,:,200)))* 12.01 * 9.0 * 60 * 60 * 24); shading flat
clim([0 5])
title('Det')

print('-dpng',[pp 'gfdl_15arcmin_2meso_test.png'])

%%
% index of water cells
[ni,nj,nt] = size(det_btm);
WID = find(~isnan(det_btm(:,:,1)));  % spatial index of water cells
NID = length(WID);                    % number of water cells

%%
for y = 2:nyrs
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
    poc= (det_btm(:,:,range));
    zm = (nmdz_100(:,:,range));
    zl = (nlgz_100(:,:,range));
    hpm = (hploss_nmdz_100(:,:,range));
    hpl = (hploss_nlgz_100(:,:,range));

    % % index of water cells
    % [ni,nj,nt] = size(tb);
    % WID = find(~isnan(tb(:,:,1)));  % spatial index of water cells
    % NID = length(WID);              % number of water cells

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

        % MZ: from gC m-2 to g(WW) m-2
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zm(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zm(j,:) = yi * 9.0;

        % LZ: from gC m-2 to g(WW) m-2
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(zl(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        Zl(j,:) = yi * 9.0;

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
    save([cpath 'Data_gfdl_mom6_cobalt2_2meso_15arcmin_daily_',num2str(YR),'.mat'], 'ESM','-v7.3');
    %     save([fpath 'gfdl_mom6_cobalt2_obsclim_15arcmin_daily_',num2str(YR),'.mat'],...
    %         'Tp','Tb','Zm','det','-v7.3');

    % TAKES A SUPER LONG TIME TO SAVE (AND FILES HUGE)
    % MAYBE RUN EACH YEAR WITH DAILY INTERP AS PART OF IT
    % INSTEAD OF SAVING

end


