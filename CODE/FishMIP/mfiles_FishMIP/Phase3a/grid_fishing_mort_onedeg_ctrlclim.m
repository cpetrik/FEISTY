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
Cdir = '/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
%Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

% Depth, lat, lon, area, grid cell with seafloor
load([Cdir 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([Cdir 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);

%  doubles
depth = double(deptho);
geolat_t = double(LAT);
geolon_t = double(LON);

[ni,nj] = size(geolat_t);

WID = GRD.ID;
NID = GRD.N;

%% need to shift lon ???
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
% take mean temp over hist time period
load([fpath 'gfdl-mom6-cobalt2_ctrlclim_temp100_onedeg_global_monthly_1961_2010.mat'],'temp_100');
load([fpath 'gfdl-mom6-cobalt2_ctrlclim_tob_onedeg_global_monthly_1961_2010.mat'],'tob');

%need to repeat 1961-1980 6x to represent 1841-1960?

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

save([fpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_mortality_all_ID_annual_tempSc.mat'],'year','WID',...
    'fmD','fmF','fmP');
save([spath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_mortality_all_ID_annual_tempSc.mat'],'year','WID',...
    'fmD','fmF','fmP');




