% Regridded fishing effort to POP grid for Additional years 2011-2015
% Check against FFmsy ending in 2010

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

%% use nominal estimate
alt1 = 'nominal';
spath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];

% added yrs for geoengineering sims
gpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/CDR_Nuclear/geoengineering_fishing_processed/';

%%
load([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_2011_2015_NOtempSc_NOmsy_',alt1,'.mat'],...
    'fmD','fmF','fmP');

fmD15 = fmD;
fmF15 = fmF;
fmP15 = fmP;
clear fmD fmF fmP

%%
load([gpath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_1948_2010_NOtempSc_NOmsy_',alt1,'.mat'],...
    'fmD','fmF','fmP')

%%
fmD2 = [fmD fmD15];
fmF2 = [fmF fmF15];
fmP2 = [fmP fmP15];

%% global mean
mD = mean(fmD2,1,'omitnan');
mF = mean(fmF2,1,'omitnan');
mP = mean(fmP2,1,'omitnan');

yrall = 1948:2015; 

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
title('forage fish enomFFmsy')
subplot(2,2,2)
plot(yrall,mP,'b')
xlim([2005 2015])
title('lg pelagics enomFFmsy')
subplot(2,2,3)
plot(yrall,mD,'k')
xlim([2005 2015])
title('demersal enomFFmsy')
print('-dpng',[gpath 'ts_nom_enomFFmsy_all.png'])

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
clim([-3 3])
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
clim([-3 3])
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
clim([-3 3])
colorbar('SouthOutside')
title('Demersals')
%print('-dpng',[gpath 'Map_diff_fmort_2011_2010_D.png'])
print('-dpng',[gpath 'Map_diff_nom_enomFFmsy_2011_2010_all.png'])

