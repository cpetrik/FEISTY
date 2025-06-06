% Visualize difference between
% CESM-WACCM Historic & CESM FOSI diatom 

clear
close all

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% CESM-WACCM
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
load([wpath 'gridspec_cesm2_cmip6_2300.mat']);
load([wpath 'Data_grid_cesm2_cmip6_2300.mat']);
load([wpath 'hist/cesm2_hist_phyc_150_monthly_1850_2014.mat']);
load([wpath 'hist/cesm2_hist_diat_150_monthly_1850_2014.mat']);

CID = GRD.ID;
CLAT = LAT;
CLON = LON;

[hi,hj]=size(CLON);

diat_waccm = double(diat_150);
phyc_waccm = double(phyc_150);
wyr = yr;

clear diat_150 phyc_150 GRD LAT LON yr

%% CESM FOSI grid
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat'],'TLAT','TLONG');
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
FID = GRD.ID;

load([cpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],'diatC_150m',...
    'diat_150m_units','time','yr','spC_150m','spC_150m_units')

diat_fosi = double(diatC_150m);
spC_fosi = double(spC_150m);
fyr=1947+(1/12):(1/12):2015;

clear diatC_150m spC_150m GRD yr

%% see if NaNs - zeros
test1 = squeeze(diat_waccm(:,:,100));
test2 = squeeze(diat_fosi(:,:,50));

figure
pcolor(test1); shading flat; colorbar

figure
pcolor(test2); shading flat; colorbar

%% convert units to gC

%FOSI zoo: nmolC cm-2 to g(WW) m-2
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
D_fosi = diat_fosi * 1e-9 * 1e4 * 12.01;
S_fosi = spC_fosi * 1e-9 * 1e4 * 12.01;
P_fosi = (D_fosi+S_fosi);

%waccm zoo: from molC m-2 to g(WW) m-2
%12.01 g C in 1 mol C
%1 g dry W in 9 g wet W (Pauly & Christiansen)
D_waccm = diat_waccm * 12.01;
P_waccm = phyc_waccm * 12.01;
S_waccm = P_waccm - D_waccm;

%% now look if same
test3 = squeeze(D_waccm(:,:,100));
test4 = squeeze(D_fosi(:,:,50));

figure
pcolor(test3); shading flat; colorbar

figure
pcolor(test4); shading flat; colorbar

%% Need to fix FOSI longitude
%tlat   [-79.2205 89.7064]
%clat   [-89.5 89.5]
%tlon   [0.0147 359.996]
%clon   [-179.5 179.5]

test = TLONG;
id=find(test>180);
test(id)=test(id)-360;
tlon = test;
tlat = TLAT;

lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

%%
figure
pcolor(tlon); shading flat; colorbar

figure
pcolor(CLON); shading flat; colorbar

figure
pcolor(glon); shading flat; colorbar

figure
pcolor(tlat); shading flat; colorbar

figure
pcolor(CLAT); shading flat; colorbar

figure
pcolor(glat); shading flat; colorbar

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

load coastlines

%% take means
[wyyr,wyid] = intersect(wyr,fyr);
[fyyr,fyid] = intersect(fyr,wyr);

mD_fosi = mean(D_fosi(:,:,fyid),3,'omitnan');
mD_waccm = mean(D_waccm(:,:,wyid),3,'omitnan');

mS_fosi = mean(S_fosi(:,:,fyid),3,'omitnan');
mS_waccm = mean(S_waccm(:,:,wyid),3,'omitnan');

mP_fosi = mean(P_fosi(:,:,fyid),3,'omitnan');
mP_waccm = mean(P_waccm(:,:,wyid),3,'omitnan');

%% Interpolate to same grid
glon=glon';
glat=glat';

fD = griddata(tlat,tlon,mD_fosi,glat,glon);
wD = griddata(CLAT,CLON,mD_waccm,glat,glon);
fS = griddata(tlat,tlon,mS_fosi,glat,glon);
wS = griddata(CLAT,CLON,mS_waccm,glat,glon);
fP = griddata(tlat,tlon,mP_fosi,glat,glon);
wP = griddata(CLAT,CLON,mP_waccm,glat,glon);

%% diffs
diffD = (wD - fD);
pdiffD = (wD-fD) ./ fD;

diffS = (wS - fS);
pdiffS = (wS-fS) ./ fS;

diffP = (wP - fP);
pdiffP = (wP-fP) ./ fP;

%%
close all

figure(1)
% WACCM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,wD)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 4]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('WACCM diat (gC)');

% FOSI
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,fD)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 4]);
colorbar
set(gcf,'renderer','painters')
title('FOSI diat');

% Diff
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,diffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 4]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI');

% PDiff
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,pdiffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI / FOSI');
%stamp('')
print('-dpng',[ppath 'WACCM_FOSI_diat_mean_global_diff_gC.png'])

%%
figure(2)
% WACCM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,wP)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('WACCM tot phy (gC)');

% FOSI
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,fP)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5]);
colorbar
set(gcf,'renderer','painters')
title('FOSI tot phy');

% Diff
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,diffP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-4 4]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI');

% PDiff
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,pdiffP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI / FOSI');
%stamp('')
print('-dpng',[ppath 'WACCM_FOSI_phy_mean_global_diff_gC.png'])


%%
figure(3)
% WACCM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,wS)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 2]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('WACCM SmP (gC)');

% FOSI
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,fS)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 2]);
colorbar
set(gcf,'renderer','painters')
title('FOSI SmP');

% Diff
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,diffS)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI');

% PDiff
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,pdiffS)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI / FOSI');
%stamp('')
print('-dpng',[ppath 'WACCM_FOSI_sp_mean_global_diff_gC.png'])
