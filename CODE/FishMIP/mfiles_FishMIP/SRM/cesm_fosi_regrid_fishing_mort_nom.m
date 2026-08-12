% Regrid fishing effort to POP grid

clear
close all

%% Map
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;
WID = GRD.ID;
NID = GRD.N;

tlme = double(lme_mask);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%% ! Shift Longitudes !
TLON = TLONG;
TLON(TLONG>180) = TLONG(TLONG>180)-360;

%% ! --> use nominal effort estimate
alt1 = 'nominal';
spath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
load([spath 'grid_mortality_all_',alt1,'.mat'])

%% ! --> subset 1948-2015
yrall = 1841:2010;
yid = find(yrall>=1948);

fmortD = fmortD(:,yid);
fmortF = fmortF(:,yid);
fmortP = fmortP(:,yid);

%% 1/2 degree
lats = unique([LatD, LatF, LatP]);
lons = unique([LonD, LonF, LonP]);

%%
nt = length(yid);
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

for t=1:nt
    clear testD testF testP

    testD = fmortD(:,t);
    D = scatteredInterpolant(LonD,LatD,testD(:));
    zGridD = D(TLON,TLAT);
    fmD(:,t) = zGridD(WID);

    testF = fmortF(:,t);
    F = scatteredInterpolant(LonF,LatF,testF(:));
    zGridF = F(TLON,TLAT);
    fmF(:,t) = zGridF(WID);

    testP = fmortP(:,t);
    P = scatteredInterpolant(LonP,LatP,testP(:));
    zGridP = P(TLON,TLAT);
    fmP(:,t) = zGridP(WID);
end

%%
clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLON,zGridF)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLON,zGridP)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLON,zGridD)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% temp scaling
load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
    'tp','tb');

[ni,nj,ft] = size(tp);

TP = reshape(tp,ni*nj,ft);
TB = reshape(tb,ni*nj,ft);

vmtp = zeros(NID,nt);
vmtb = zeros(NID,nt);

for t=1:nt
    clear testB testP

    testB = TB(:,t);
    vmtb(:,t) = testB(WID);

    testP = TP(:,t);
    vmtp(:,t) = testP(WID);
end

% vmtp = mtp(WID);
% vmtb = mtb(WID);

%% save before scaling
year = 1948:2010;

gpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/CDR_Nuclear/geoengineering_fishing_processed/';

save([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_1948_2010_NOtempSc_NOmsy_',alt1,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

%% scale with Fmsy and temp ! --> find temp means from FOSI !
%fm = F/Fmsy, need to mult by Fmsy ~= 0.3
%tsc = (exp(0.063*(temp-10.0)));

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
year = 1948:2010;

gpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/CDR_Nuclear/geoengineering_fishing_processed/';

save([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_1948_2010_tempSc_',alt1,'.mat'],'year','WID',...
    'fmD','fmF','fmP');



