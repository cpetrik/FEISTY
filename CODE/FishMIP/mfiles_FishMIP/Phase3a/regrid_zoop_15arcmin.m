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

%% Are they the same orientation? - No

lat2 = flipud(lat);
lon2 = flipud(lon);
lat3 = lat2';
lon3 = lon2';

figure(1)
pcolor(lat2'); shading flat;
colorbar
title('lat')

figure(2)
pcolor(LAT); shading flat;
colorbar
title('LAT')

figure(3)
pcolor(lon2'); shading flat;
colorbar
title('lon')

figure(4)
pcolor(LON); shading flat;
colorbar
title('LON')

figure(5)
pcolor(LON,LAT,deptho); shading flat
colorbar
title('depth')

%% shift lon
lon4 = lon3 + 119.75;

figure(6)
pcolor(lon4); shading flat;
colorbar
title('lon4')

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

%% Test on depth
% F = scatteredInterpolant(x,y,v); % fn of old grid
% vq = F(xq,yq);                   % interp on new grid
% F = griddedInterpolant(x,y,v);   % fn of old grid
% vq = F(xq,yq);                   % interp on new grid

F = griddedInterpolant(fliplr(LON),fliplr(LAT),fliplr(deptho)); % fn of old grid
dep2 = F(fliplr(lon4),fliplr(lat3));                   % interp on new grid

%%
close all
figure(7)
pcolor(LON,LAT,deptho); shading flat
colorbar
title('depth')

figure(8)
pcolor(lon4,lat3,fliplr(dep2)); shading flat
colorbar
title('depth2')

%% Test going in opposite direction - THIS THE WAY
dep2 = fliplr(dep2);
dep3 = log10(abs(dep2) +1);
D = griddedInterpolant(fliplr(lon4),fliplr(lat3),fliplr(dep3)); % fn of old grid
dep4 = D(fliplr(LON),fliplr(LAT));                   % interp on new grid

%%
close all
figure(9)
pcolor(LON,LAT,fliplr(dep4)); shading flat
colorbar
title('dep4')

figure(10)
pcolor(lon4,lat3,dep3); shading flat
colorbar
title('dep3')


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


