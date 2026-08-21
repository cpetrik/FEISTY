% Put fishing mortality onto grid

clear 
close all

%% 1993-2010
% Hopefully get 2011-2015

% Predators use v1
alt1 = 'grid_mortality_guilds_v1';

spath1 = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt1,'/'];
load([spath1 'grid_mortality_all_v1.mat'],'fmortD','fmortP','LatD','LatP',...
    'LonD','LonP','yrD','yrP')

% Forage use v3
alt3 = 'grid_mortality_guilds_v3';
spath3 = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt3,'/'];
load([spath3 'grid_mortality_all_v3.mat'],'fmortF','LatF','LonF','yearF')

spath32 = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_impl/grid_mortality_guilds_v32/'];

%% 1841-2010, subset 1961-2010
yrall = 1841:2010;
yid = find(yrall>=1993);

fmortD = fmortD(:,yid);
fmortF = fmortF(:,yid);
fmortP = fmortP(:,yid);

%% 1/2 degree
lats = unique([LatD, LatF, LatP]);
lons = unique([LonD, LonF, LonP]);

%% NEP grid info
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';

% Depth, lat, lon, area, grid cell with seafloor
load([cpath 'nep_raw_ocean_static_gridspec.mat'],'geolon','geolat');
load([cpath 'Data_grid_mom6_nep10.mat'],'GRD');

geolon = double(geolon);
geolat = double(geolat);
geolon(geolon>1e20) = nan;
geolat(geolat>1e20) = nan;

[ni,nj] = size(geolon);

WID = GRD.ID;
NID = GRD.NID;

%%
figure
pcolor(geolon'); shading flat; colorbar

figure
pcolor(geolat'); shading flat; colorbar

%% make W longitudes neg
sid = find(geolon>180);
geolon(sid) = geolon(sid) - 360;

%%
nt = length(yid);
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

for t=1:nt
    clear testD testF testP
    
    testD = griddata(LonD,LatD,fmortD(:,t),geolon,geolat);
    fmD(:,t) = testD(WID);
    
    testF = griddata(LonF,LatF,fmortF(:,t),geolon,geolat);
    fmF(:,t) = testF(WID);
    
    testP = griddata(LonP,LatP,fmortP(:,t),geolon,geolat);
    fmP(:,t) = testP(WID);
end

%%
%NE Pac
plotminlat=10; 
plotmaxlat=85;
plotminlon=156;
plotmaxlon=-104;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

figure
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,testF)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,testP)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','off',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,testD)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% temp scaling
load([cpath 'nep.full.hcast.monthly.raw.r20250912.199301-202506_annual_means.mat'],...
    'tp_amean','tb_amean','yrs')

vmtp = nan*ones(NID,length(yrs));
vmtb = nan*ones(NID,length(yrs));

for y=1:length(yrs)
    tp = tp_amean(:,:,t);
    tb = tb_amean(:,:,t);
    vmtp(:,y) = tp(WID);
    vmtb(:,y) = tb(WID);
end

%% scale with Fmsy and temp
% fm = F/Fmsy, need to mult by Fmsy ~= 0.3
%tsc = (exp(0.063*(temp-10.0));

Fyr = yrall(yid);
[tyr,tid] = intersect(Fyr,yrs);

fmF = 0.3 * fmF .* (exp(0.063*(vmtp(:,tid)-10.0)));
fmP = 0.3 * fmP .* (exp(0.063*(vmtp(:,tid)-10.0)));
fmD = 0.3 * fmD .* (exp(0.063*(vmtb(:,tid)-10.0)));

%%
fmD(isnan(fmD)) = 0.0;
fmF(isnan(fmF)) = 0.0;
fmP(isnan(fmP)) = 0.0;

fmD(fmD<0) = 0.0;
fmF(fmF<0) = 0.0;
fmP(fmP<0) = 0.0;

%% save 
year = 1993:2010;

save([cpath 'nep.full.hcast.monthly.raw.r20250912_fmort_annual_1993_2010_tempSc_v32.mat'],'year','WID',...
    'fmD','fmF','fmP');
save([spath32 'nep.full.hcast.monthly.raw.r20250912_fmort_annual_1993_2010_tempSc_v32.mat'],'year','WID',...
    'fmD','fmF','fmP');




