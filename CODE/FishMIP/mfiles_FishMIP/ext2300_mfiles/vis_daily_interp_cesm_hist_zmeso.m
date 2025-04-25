% View CESM mesozoo
% Hist 1850
% Using 100% diat, 43% SP fractions

clear 
close all

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% Units
%poc flux: mol C m-2 s-1
%zoo: mol C m-2
%tp: degC
%tb: degC

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
load([fpath 'cesm2_hist_zmeso_150_monthly_1850_2014.mat']);
load([fpath 'Data_cesm_hist_daily_1850.mat']);

cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
load([cpath 'gridspec_cesm2_cmip6_2300.mat']);
load([cpath 'Data_grid_cesm2_cmip6_2300.mat']);

CESM = ESM;
CID = GRD.ID;
CLAT = LAT;
CLON = LON;

clear ESM GRD LAT LON

%%
ipath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/hist/';
load([ipath 'Data_ipsl_hist_daily_1850.mat']);
ipath2 = ['/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/'];

load([ipath2 'gridspec_ipsl_cmip6_2300.mat']);
load([ipath2 'Data_grid_ipsl_cmip6_2300.mat']);

IPSL = ESM;
IID = GRD.ID;
ILAT = LAT;
ILON = LON;

clear ESM GRD LAT LON

%%
upath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/hist/';
load([upath 'Data_ukesm_hist_daily_1850.mat']);
upath2 = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/';
  
    load([upath2 'gridspec_ukesm_cmip6_2300.mat']);
    load([upath2 'Data_grid_ukesm_cmip6_2300.mat']);

UK = ESM;
UID = GRD.ID;
ULAT = LAT;
ULON = LON;

clear ESM GRD LAT LON

%% Map data

[ni,nj]=size(CLON);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines

%%
Zm = mean(zmeso_150(:,:,1:12),3,'omitnan') * 12.01 * 9.0;

%%
pcolor(log10(Zm))
shading flat
cmocean('tempo')
clim([0 2])
colorbar
title([mod ' Zmeso43 (gWW/m2)'])

%%
C_Zm = mean(CESM.Zm,2,'omitnan');
I_Zm = mean(IPSL.Zm,2,'omitnan');
U_Zm = mean(UK.Zm,2,'omitnan');

Czm = NaN*ones(ni,nj);
Izm = NaN*ones(ni,nj);
Uzm = NaN*ones(ni,nj);

Czm(CID) = C_Zm;
Izm(IID) = I_Zm;
Uzm(UID) = U_Zm;

ebins = -1:0.25:2;

%%

f1 = figure('Units','inches','Position',[1 3 5 5]);

subplot('Position',[0.05 0.68 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,log10(Czm))
cmocean('tempo')
clim([0 2])
colorbar
title('CESM Zmeso (gWW/m2)')
%patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.68 0.425 0.275])
histogram(log10(Czm(:)),ebins)

subplot('Position',[0.05 0.37 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(ILAT,ILON,log10(Izm))
cmocean('tempo')
clim([0 2])
colorbar
title('IPSL Zmeso (gWW/m2)')
%h3=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.37 0.425 0.275])
histogram(log10(Izm(:)),ebins)

subplot('Position',[0.05 0.06 0.45 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(ULAT,ULON,log10(Uzm))
cmocean('tempo')
clim([0 2])
colorbar
title('UKESM Zmeso (gWW/m2)')
%h5=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.06 0.425 0.275])
histogram(log10(Uzm(:)),ebins)

print('-dpng',[ppath 'map_ESMs_hist_zmeso_1850_Lfrac30.png'])

%% Plot t.s. with other ESMs
Uhist_Zm = mean(UK.Zm,'omitnan');
Ihist_Zm = mean(IPSL.Zm,'omitnan');
Chist_Zm = mean(CESM.Zm,'omitnan');

%%

figure
plot(1:365,Uhist_Zm,'b','LineWidth',2); hold on
plot(1:365,Ihist_Zm,'r','LineWidth',2); hold on
plot(1:365,Chist_Zm,'color',[0 0.75 0.5],'LineWidth',2);
title('Mesozoo')
legend('UKESM','IPSL','CESM43')
legend('location','northwest')
print('-dpng',[ppath 'ts_ESMs_hist_zmeso_1850_Lfrac30.png'])



