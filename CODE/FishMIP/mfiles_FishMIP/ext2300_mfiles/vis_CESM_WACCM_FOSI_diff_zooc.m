% Visualize difference between
% CESM-WACCM Historic & CESM FOSI zooc

clear
close all

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% CESM-WACCM
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
load([wpath 'gridspec_cesm2_cmip6_2300.mat']);
load([wpath 'Data_grid_cesm2_cmip6_2300.mat']);
load([wpath 'hist/cesm2_hist_zooc_150_monthly_1850_2014.mat']);

CID = GRD.ID;
CLAT = LAT;
CLON = LON;

[hi,hj]=size(CLON);

z_waccm = double(zooc_150);
wyr = yr;

clear zooc_150 GRD LAT LON yr

%% CESM FOSI grid
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat'],'TLAT','TLONG');
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
FID = GRD.ID;

load([cpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],'zooC_150m',...
    'zooC_150m_units','time','yr')

z_fosi = double(zooC_150m);
fyr=1947+(1/12):(1/12):2015;

clear zooC_150m GRD yr

%% see if NaNs - no
test1 = squeeze(z_waccm(:,:,100));
test2 = squeeze(z_fosi(:,:,50));

figure
pcolor(test1); shading flat; colorbar

figure
pcolor(test2); shading flat; colorbar

%% convert units to gWW

%FOSI zoo: nmolC cm-2 to g(WW) m-2
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
Z_fosi = z_fosi * 1e-9 * 1e4 * 12.01 * 9.0;

%waccm zoo: from molC m-2 to g(WW) m-2
%12.01 g C in 1 mol C
%1 g dry W in 9 g wet W (Pauly & Christiansen)
Z_waccm = z_waccm * 12.01 * 9.0;

%% now look if same
test3 = squeeze(Z_waccm(:,:,100));
test4 = squeeze(Z_fosi(:,:,50));

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

mZ_fosi = mean(Z_fosi(:,:,fyid),3,'omitnan');
mZ_waccm = mean(Z_waccm(:,:,wyid),3,'omitnan');

%% Interpolate to same grid
glon=glon';
glat=glat';

fZ = griddata(tlat,tlon,mZ_fosi,glat,glon);

wZ = griddata(CLAT,CLON,mZ_waccm,glat,glon);

%% diffs
diffZ = (wZ - fZ);
pdiffZ = (wZ-fZ) ./ fZ;

wZ_gC = wZ/9;
fZ_gC = fZ/9;
CdiffZ = (wZ_gC - fZ_gC);
CpdiffZ = (wZ_gC-fZ_gC) ./ fZ_gC;

%%
close all

figure(10)
% WACCM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,wZ)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 15]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('WACCM zooc (gWW)');

% FOSI
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,fZ)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 15]);
colorbar
set(gcf,'renderer','painters')
title('FOSI zooc');

% Diff
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,diffZ)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-8 8]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI');

% PDiff
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,pdiffZ)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI / FOSI');
%stamp('')
print('-dpng',[ppath 'WACCM_FOSI_zooc_mean_global_diff_gWW.png'])

%%
figure(11)
% WACCM
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,wZ_gC)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 2]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('WACCM zooc (gC)');

% FOSI
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,fZ_gC)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 2]);
colorbar
set(gcf,'renderer','painters')
title('FOSI zooc');

% Diff
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,CdiffZ)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.8 0.8]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI');

% PDiff
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,CpdiffZ)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('WACCM - FOSI / FOSI');
%stamp('')
print('-dpng',[ppath 'WACCM_FOSI_zooc_mean_global_diff_gC.png'])