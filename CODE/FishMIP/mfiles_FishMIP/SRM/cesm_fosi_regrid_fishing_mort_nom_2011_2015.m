% Regrid fishing effort to POP grid
% Additional years 2011-2015

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

%% ! --> use nominal estimate
alt1 = 'nominal';
spath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
% load([spath 'grid_mortality_all_',alt1,'.mat'])

% added yrs for geoengineering sims
gpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/CDR_Nuclear/geoengineering_fishing_processed/';
load([gpath 'all_enom_FFmsy_2011-2015.mat'])

%% 1/2 degree
lats = unique([DLat; FLat; PLat]);
lons = unique([DLon; FLon; PLon]);
yid = unique([FYear;PYear;DYear]);

%%
nt = length(yid);
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

for t=1:nt
    clear testD testF testP

    did = find(DYear==yid(t));
    D = scatteredInterpolant(DLon(did),DLat(did),demFFmsy(did));
    zGridD = D(TLON,TLAT);
    fmD(:,t) = zGridD(WID);

    fid = find(FYear==yid(t));
    F = scatteredInterpolant(FLon(fid),FLat(fid),forageFFmsy(fid));
    zGridF = F(TLON,TLAT);
    fmF(:,t) = zGridF(WID);

    pid = find(PYear==yid(t));
    P = scatteredInterpolant(PLon(pid),PLat(pid),lgpelFFmsy(pid));
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
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLON,zGridP)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLON,zGridD)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% temp scaling
load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
    'tp','tb');

yrall = 1948:2015; 

[ni,nj,ft] = size(tp);

TP = reshape(tp,ni*nj,ft);
TB = reshape(tb,ni*nj,ft);

vmtp = zeros(NID,nt);
vmtb = zeros(NID,nt);

for t=1:nt
    clear testB testP

    yr = find(yrall==yid(t));

    testB = TB(:,yr);
    vmtb(:,t) = testB(WID);

    testP = TP(:,yr);
    vmtp(:,t) = testP(WID);
end

% vmtp = mtp(WID);
% vmtb = mtb(WID);

%% save before temp scaling
year = 2011:2015;

save([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_2011_2015_NOtempSc_NOmsy_',alt1,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

%% scale with Fmsy and temp 
%fm = F/Fmsy, need to mult by Fmsy ~= 0.3
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
year = 2011:2015;

save([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_2011_2015_tempSc_',alt1,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

%% Check for jumps between 2010 and 2011
fmD15 = fmD;
fmF15 = fmF;
fmP15 = fmP;
clear fmD fmF fmP

load([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_1948_2010_tempSc_',alt1,'.mat'])

%%
fmD2 = [fmD fmD15];
fmF2 = [fmF fmF15];
fmP2 = [fmP fmP15];

%% global mean
mD = mean(fmD2,1,'omitnan');
mF = mean(fmF2,1,'omitnan');
mP = mean(fmP2,1,'omitnan');

figure(1)
subplot(2,2,1)
plot(yrall,mF,'r')
subplot(2,2,2)
plot(yrall,mP,'b')
subplot(2,2,3)
plot(yrall,mD,'k')

figure(2)
subplot(2,2,1)
plot(yrall,mF,'r')
xlim([2005 2015])
title('forage fish fmort rate')
subplot(2,2,2)
plot(yrall,mP,'b')
xlim([2005 2015])
title('lg pelagics fmort rate')
subplot(2,2,3)
plot(yrall,mD,'k')
xlim([2005 2015])
title('demersal fmort rate')
print('-dpng',[gpath 'ts_fishing_mort_nom_convert_enomFFmsy_all.png'])

%% global diff 2011 from 2010
dF = fmF15(:,1) - fmF(:,end);
dP = fmP15(:,1) - fmP(:,end);
dD = fmD15(:,1) - fmD(:,end);

gF = nan(ni,nj);
gP = nan(ni,nj);
gD = nan(ni,nj);

gF(WID) = dF;
gP(WID) = dP;
gD(WID) = dD;

figure(4)
subplot('Position',[0.05 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,gF)
hold on
cmocean('balance')
clim([-1 1])
colorbar('SouthOutside')
title('Forage Fish')
%print('-dpng',[gpath 'Map_diff_fmort_2011_2010_F.png'])

subplot('Position',[0.5 0.5 0.45 0.45])
%figure(5)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,gP)
hold on
cmocean('balance')
clim([-1 1])
colorbar('SouthOutside')
title('Large Pelagics')
%print('-dpng',[gpath 'Map_diff_fmort_2011_2010_P.png'])

%figure(6)
subplot('Position',[0.05 0.05 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,gD)
hold on
cmocean('balance')
clim([-1 1])
colorbar('SouthOutside')
title('Demersals')
%print('-dpng',[gpath 'Map_diff_fmort_2011_2010_D.png'])
print('-dpng',[gpath 'Map_diff_fmort_nom_2011_2010_all.png'])



