% Make mat files of interpolated time series from GFDL
% CORE-forced runs 1988-2007
% velocity means top 100m

clear
close all

gpath = '/project/Feisty/GCM_DATA/CORE-forced/';

%%
load([gpath 'ocean_cobalt_grid.mat'],'geolon_t','geolat_t');
load([gpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);

%% 
load([gpath 'vel100_feb152013_run25_ocean.198801-200712.mat']);

%%
mos = length(TIME);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;

yrs = floor(yr(1)):yr(end);

Tdays=1:365;

%% land mask

lmask = nan*ones(ni,nj);
lmask(GRD.ID) = GRD.lmask;
mask = repmat(lmask,1,1,mos);

u_100 = u_100 .* mask;
v_100 = v_100 .* mask;

%% test that all same orientation
load([gpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');

tp = nan*ones(ni,nj);
tp(GRD.ID) = COBALT.Tp(:,1);

%%
close all

% index of water cells
WID = GRD.ID;
NID = length(WID);

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

    u = (u_100(:,:,range));
    v = (v_100(:,:,range));

    % setup FEISTY data files
    U100  = nan*zeros(NID,365);
    V100  = nan*zeros(NID,365);

    % interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell

        % pelagic temperature (in Celcius)
        X = squeeze(u(m,n,:));
        xi = interp1(Time, X, 1:365,'linear','extrap');
        U100(j,:) = xi;

        % bottom temperature (in Celcius)
        Y = squeeze(v(m,n,:));
        yi = interp1(Time, Y, 1:365,'linear','extrap');
        V100(j,:) = yi;

    end

    ESM.U  = U100;
    ESM.V  = V100;

    % save
    save([gpath 'Vel100_esm2m_core_daily_',num2str(YR),'.mat'],'ESM','-v7.3');
    
end


    
