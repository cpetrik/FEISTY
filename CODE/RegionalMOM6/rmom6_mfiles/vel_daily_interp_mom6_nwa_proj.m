% Make mat files of interpolated time series from GFDL
% Reanalysis-forced runs 1993-2019
% MOM6-NWA12 velocity means top 100m

clear
close all

fpath='/project/Feisty/GCM_Data/MOM6-NWA12/';

%% 
load([fpath 'u_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],'u_100');
load([fpath 'v_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat']);

%%
u_100(u_100 > 1.0e15) = nan;
v_100(v_100 > 1.0e15) = nan;

%%
load([fpath 'nwa_raw_ocean_static_gridspec.mat'],'geolon','geolat');

geolat = double(geolat);
geolon = double(geolon);

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

Tdays=1:365;

%% need to regrid velocities to tracer grid
[ni,nj] = size(geolon);

% [~,~,nt] = size(u_100);
% 
% ut_100 = nan*ones(ni,nj,nt);
% vt_100 = nan*ones(ni,nj,nt);
% 
% for t=1:nt
%     u = u_100(:,:,t);
%     v = v_100(:,:,t);
% 
%     ut = griddata(geolat_u,geolon_u,u,geolat,geolon);
%     vt = griddata(geolat_v,geolon_v,v,geolat,geolon);
% 
%     ut_100(:,:,t) = ut;
%     vt_100(:,:,t) = vt;
% end

% Introduces NaNs, just remove one row or col instead
ut_100 = u_100(1:ni,:,:);
vt_100 = v_100(:,2:end,:);

%%
clear u_100 v_100 u v ut vt

close all

%%
load([fpath 'Data_grid_mom6_nwa12.mat'], 'GRD');

% index of water cells
UVID = GRD.UVID;
NID = length(UVID);

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

    u = (ut_100(:,:,range));
    v = (vt_100(:,:,range));

    % setup FEISTY data files
    U100  = nan*zeros(NID,365);
    V100  = nan*zeros(NID,365);

    % interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],UVID(j)); % spatial index of water cell

        % pelagic temperature (in Celcius)
        X = squeeze(u(m,n,:));
        xi = interp1(Time, X, 1:365,'linear','extrap');
        U100(j,:) = xi;

        % bottom temperature (in Celcius)
        Y = squeeze(v(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        V100(j,:) = yi;

    end

    MOM.U  = U100;
    MOM.V  = V100;

    % save
    save([fpath 'Data_mom6_nwa12_daily_',num2str(YR),'.mat'],'MOM','-v7.3');
    
end


    
