% Put fishing mortality onto grid

clear all 
close all

%% 1961-2010
spath = '/Volumes/MIP/Fish-MIP/Phase3/fishing/grid_mortality_guilds/';
fpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/grid_mortality_guilds/';
load([fpath 'grid_mortality_all.mat'])

%% 1/2 degree
lats = unique([LatD, LatF, LatP]);
lons = unique([LonD, LonF, LonP]);

%% test with CORE-forced COBALT
Cdir = '/Volumes/MIP/GCM_DATA/CORE-forced/';
%Cdir = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
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
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

for t=1:nt
    testD = griddata(LonD,LatD,fmortD(:,t),geolon_t,geolat_t);
    fmD(:,t) = testD(WID);
    
    testF = griddata(LonF,LatF,fmortF(:,t),geolon_t,geolat_t);
    fmF(:,t) = testF(WID);
    
    testP = griddata(LonP,LatP,fmortP(:,t),geolon_t,geolat_t);
    fmP(:,t) = testP(WID);
    
    clear testD testF testP
end

%% temp scaling
load([Cdir 'CORE_interann_mean_forcings_anom.mat'],'tp','tb');

mtp = nanmean(tp,3);
mtb = nanmean(tb,3);

vmtp = mtp(WID);
vmtb = mtb(WID);

%%
%tsc = (exp(0.063*(temp-10.0));

fmF = 0.3 * fmF .* (exp(0.063*(vmtp-10.0)));
fmP = 0.3 * fmP .* (exp(0.063*(vmtp-10.0)));
fmD = 0.3 * fmD .* (exp(0.063*(vmtb-10.0)));

%shelf
sid = find(GRD.Z<200);
%deep
did = find(GRD.Z>=200);

%fmD(did,:) = fmD(did,:) .* vmtb(did);

%%
fmD(isnan(fmD)) = 0.0;
fmF(isnan(fmF)) = 0.0;
fmP(isnan(fmP)) = 0.0;

fmD(fmD<0) = 0.0;
fmF(fmF<0) = 0.0;
fmP(fmP<0) = 0.0;

%% save 
year = 1961:2010;

save([fpath 'CORE_mortality_all_ID_annual_tempSc.mat'],'year','WID',...
    'fmD','fmF','fmP');
save([spath 'CORE_mortality_all_ID_annual_tempSc.mat'],'year','WID',...
    'fmD','fmF','fmP');




