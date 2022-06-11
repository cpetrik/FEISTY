% Put fishing mortality onto grid

clear all 
close all

%% 1961-2010
fpath = '/Volumes/MIP/Fish-MIP/Phase3/fishing/grid_mortality_guilds/';
load([fpath 'grid_mortality_all.mat'])

%% 1/2 degree
lats = unique([LatD, LatF, LatP]);
lons = unique([LonD, LonF, LonP]);

%% test with CORE-forced COBALT
Cdir = '/Volumes/MIP/GCM_DATA/CORE-forced/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

% Depth, lat, lon, area, grid cell with seafloor
load([Cdir 'ocean_cobalt_grid.mat'])

%  doubles
kmt = double(kmt);
depth = double(ht);
geolat_t = double(geolat_t);
geolon_t = double(geolon_t);
area = double(area_t);

[ni,nj] = size(geolat_t);

load([Cdir 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
WID = GRD.ID;
NID = GRD.N;

%% need to shift lon
%fix lon shift
id=find(geolon_t(:)<-180);
geolon_t(id)=geolon_t(id)+360;

%%
nt = length(yearD);
fmD = nan*ones(NID,nt);
fmF = nan*ones(NID,nt);
fmP = nan*ones(NID,nt);

for t=1:nt
    testD = griddata(LonD,LatD,fmortD(:,t),geolon_t,geolat_t);
    fmD(:,t) = testD(WID);
    
    testF = griddata(LonF,LatF,fmortF(:,t),geolon_t,geolat_t);
    testF(land) = nan;
    fmF(:,t) = testF(WID);
    
    testP = griddata(LonP,LatP,fmortP(:,t),geolon_t,geolat_t);
    testP(land) = nan;
    fmP(:,t) = testP(WID);
    
    clear testD testF testP
end

%% save 
year = 1961:2010;

save([fpath 'CORE_mortality_all_ID_annual.mat'],'year','WID',...
    'fmD','fmF','fmP');




