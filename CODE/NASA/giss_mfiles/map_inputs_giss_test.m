% CESM FOSI output 

clear 
close all

%% Paths

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/VolMIP/';

load([fpath 'giss_10yr_FEISTY_forcing.mat']);

load([fpath 'giss_10yr_gridspec.mat']);

%%
[LAT,LON] = meshgrid(lat,lon);

%% means
cTp = mean(ptemp_upper200m_avg,3,'omitnan');
cTb = mean(ptemp_bottom,3,'omitnan');
cD = mean(poc_mgcm3_bottom,3,'omitnan');
cZ = mean(biomass_Zmeso_intrg200m,3,'omitnan');

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,cTp)
cmocean('thermal')
clim([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,cTb)
cmocean('thermal')
clim([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cZ))
cmocean('tempo')
clim([0 4])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zmeso (kgC m^-^2)')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(cD))
cmocean('tempo')
clim([0 2])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det (mgC cm^-^3)')
print('-dpng',[pp 'Map_giss_10yr_mean_forcings.png'])

%% convert units
% wgs84 = wgs84Ellipsoid("km");
% grid_area = areaint(LAT,LON,wgs84);

% Zmeso kgC -> gWW
% 1e3 gC/kgC * 9 gWW/gC

% Det mgC/cm3 -> gWW/m2
% 1e-3 gC/mgC * 9 gWW/gC
% cm3 * cm_gridcell_thickness -> cm2
% cm2 * 1e-4m2/cm


