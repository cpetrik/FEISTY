% regrid 1/4 files form Xiao to regular 1/4 grid
% to match up with ISIMIP vars for tp, tb, det, etc. 
% Xiao: 1440x1080, ISIMIP: 1440x720

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';


load([fpath 'gfdl-mom6-cobalt2_15arcmin_mdz_lgz_100m_global_monthly_1961_2010.mat'])
load([fpath 'gfdl-mom6_cobalt2_15arcmin_HPloss_mdz_lgz_100m_month_1961_2010.mat'])

load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
load([cpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])
%[cpath 'gfdl-mom6-cobalt2_obsclim_zmeso100_15arcmin_global_monthly_1961_2010.mat'];
%load([cpath 'Data_gfdl_mom6_cobalt2_obsclim_15arcmin_daily_2010.mat'])

%% Are they the same orientation? - No: flipud & transpose

lat2 = flipud(lat);
lon2 = flipud(lon);
lat3 = lat2';
lon3 = lon2';

% shift lon
lat4 = lat3;
lon4 = lon3;
lon4(lon3<-180) = lon3(lon3<-180) + 360; %lon4 = lon3 + 119.75;

%%
figure(1)
pcolor(lat4); shading flat;
colorbar
title('lat4')

figure(2)
pcolor(LAT); shading flat;
colorbar
title('LAT')

figure(3)
pcolor(lon4); shading flat;
colorbar
title('lon4')

figure(4)
pcolor(LON); shading flat;
colorbar
title('LON')

figure(5)
pcolor(lon4); shading flat;
colorbar
title('lon4')

figure(6)
pcolor(deptho); shading flat
colorbar
title('depth')

%% THIS IS THE WAY
% D = griddedInterpolant(fliplr(lon4),fliplr(lat3),fliplr(dep3)); % fn of old grid
% dep4 = D(fliplr(LON),fliplr(LAT));                   % interp on new grid

%%
t=1; 
zm = nmdz_100(:,:,t);
zm2 = fliplr(zm);    % flippingLR matches orientation of lat4, LAT & LON Deptho
% zm2 = flipud(zm);
% zm3 = zm2';

figure
pcolor(zm); shading flat
colorbar
title('zm')

figure
pcolor(zm2); shading flat
colorbar
title('zm2')

%% Sort lat and lon in ascending order
% flipLR both to get lat ascending
lat5 = fliplr(lat4);
lon5 = fliplr(lon4);
zm5 = zm;

% need to move rows of lon to get ascending
%neg lon starts at row 484
lon6 = lon5(484:end,:);
lon6(958:1440,:) = lon5(1:483,:);

lat6 = lat5(484:end,:);
lat6(958:1440,:) = lat5(1:483,:);

zm6 = zm5(484:end,:);
zm6(958:1440,:) = zm5(1:483,:);

%%
figure
pcolor(lat6); shading flat
colorbar
title('lat6')

figure
pcolor(lon6); shading flat
colorbar
title('lon6')

figure
pcolor(zm6); shading flat
colorbar
title('zm6')

%% Interp
%Error using griddedInterpolant
%Sample points must be sorted in ascending order.

F1 = griddedInterpolant(lon6,lat6,zm6);
test1 = F1(fliplr(LON),fliplr(LAT));

%%
close all
figure(9)
pcolor(LON,LAT,fliplr(test1)); shading flat
colorbar
title('test1')

figure(10)
pcolor(LON,LAT,deptho); shading flat
colorbar
title('deptho')

lidZ = find(isnan(fliplr(test1)));
lidD = find(isnan(fliplr(deptho)));

%Successfully regridded, but 1/4 grids do not align !!!!!!!!!

%%
[ni,nj] = size(lon4); %lat,lon 1080x1440
[zi,zj] = size(LON); %LAT,LON 1440x720

WID = GRD.ID;
NID = GRD.N;

%%
nt = length(yr);
MZ = zeros(zi,zj,nt);
LZ = MZ;
hpMZ = MZ;
hpLZ = MZ;

%%

for t=1:nt
    clear testD testF testP
    
    %test1 = griddata(LonD,LatD,nmdz_100(:,:,t),LON,LAT);
    F1 = griddedInterpolant(LonD,LatD,nmdz_100(:,:,t));
    test1 = F1(LON,LAT);
    MZ(:,:,t) = test1(WID);
    
    %test2 = griddata(LonF,LatF,nlgz_100(:,:,t),LON,LAT);
    LZ(:,:,t) = test2(WID);
    
    %test3 = griddata(LonP,LatP,hploss_nmdz_100(:,:,t),LON,LAT);
    hpMZ(:,:,t) = test3(WID);

    %test4 = griddata(LonP,LatP,hploss_nlgz_100(:,:,t),LON,LAT);
    hpLZ(:,:,t) = test4(WID);
end


