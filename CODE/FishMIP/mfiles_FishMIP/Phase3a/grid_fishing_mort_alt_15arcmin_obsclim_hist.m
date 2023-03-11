% Put fishing mortality onto grid

clear 
close all

%% 1961-2010
alt1 = 'grid_mortality_guilds_v3'; %grid_mortality_guilds_v2, pristine_grid_mortality_guilds_v1
alt2 = '_v3'; %_v2, _v1_pristine
spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
fpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt1,'/'];
load([fpath 'grid_mortality_all',alt2,'.mat'])

%% 1841-2010, subset 1961-2010
yrall = 1841:2010;
yid = find(yrall>=1961);

fmortD = fmortD(:,yid);
fmortF = fmortF(:,yid);
fmortP = fmortP(:,yid);

%% 1/2 degree
lats = unique([LatD, LatF, LatP]);
lons = unique([LonD, LonF, LonP]);

%% obsclim one degree
Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

% Depth, lat, lon, area, grid cell with seafloor
load([Cdir 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat']);
load([Cdir 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat']);

[ni,nj] = size(LON);

WID = GRD.ID;
NID = GRD.N;

%%
nt = length(yid);
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

for t=1:nt
    clear testD testF testP
    
    testD = griddata(LonD,LatD,fmortD(:,t),LON,LAT);
    fmD(:,t) = testD(WID);
    
    testF = griddata(LonF,LatF,fmortF(:,t),LON,LAT);
    fmF(:,t) = testF(WID);
    
    testP = griddata(LonP,LatP,fmortP(:,t),LON,LAT);
    fmP(:,t) = testP(WID);
end

%%
clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testF)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testP)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testD)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% temp scaling
tpath ='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([tpath 'gfdl-mom6-cobalt2_obsclim_mtemp_15arcmin_global_1961_2010.mat'])
vmtp = mtp(WID);
vmtb = mtb(WID);

%% scale with Fmsy and temp
% fm = F/Fmsy, need to mult by Fmsy ~= 0.3
%tsc = (exp(0.063*(temp-10.0));

fmF = 0.3 * fmF .* (exp(0.063*(vmtp-10.0)));
fmP = 0.3 * fmP .* (exp(0.063*(vmtp-10.0)));
fmD = 0.3 * fmD .* (exp(0.063*(vmtb-10.0)));

%%
fmD(isnan(fmD)) = 0.0;
fmF(isnan(fmF)) = 0.0;
fmP(isnan(fmP)) = 0.0;

fmD(fmD<0) = 0.0;
fmF(fmF<0) = 0.0;
fmP(fmP<0) = 0.0;

%% save 
year = 1961:2010;

% save([tpath 'gfdl-mom6-cobalt2_obsclim_15arcmin_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
%     'fmD','fmF','fmP');
save([fpath 'gfdl-mom6-cobalt2_obsclim_15arcmin_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
    'fmD','fmF','fmP');
save([spath 'gfdl-mom6-cobalt2_obsclim_15arcmin_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
    'fmD','fmF','fmP');




