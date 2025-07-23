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
load([wpath 'hist/cesm2_hist_zmeso_150_diatfraconly_monthly_1850_2014.mat'],'zmeso_150');

CID = GRD.ID;
CLAT = LAT;
CLON = LON;

[ni,nj,nt] = size(zmeso_150);

clear LAT LON 

%% GLMM of mesozoo from Heneghan
gpath = '/Users/cpetrik/Dropbox/Fish-MIP/CMIP6/driver_analysis/zmeso_matlab/';

load([gpath 'cmip6_hist_space_means_50yr_zmeso200_glmm100_same_orientation.mat'],...
    'ozmo_all','lat_g','lon_g');
 
%obsglm mg C/m2 
%50-year period from 1965 to 2014

%%
figure
pcolor(squeeze(zmeso_150(:,:,end)))
shading flat
title('zmeso')

figure
pcolor(ozmo_all)
shading flat
title('obsGLMM')

figure
pcolor(lat_g)
shading flat
title('lat_g')

figure
pcolor(CLAT)
shading flat
title('Clat')

%% flip L-R GLMM
close all

glmm = fliplr(ozmo_all);
GLAT = fliplr(lat_g);
GLON = fliplr(lon_g);

%% convert units to molC m-2

%obsglm: mg C/m2 
% 1 gC in 1e-3 mgC
% 1 molC in 12.01 g C
glmm = glmm * 1e-3 * (1/12.0);

%% 50-yr mean of CESM
y50 = find(yr>1965);
mzooc = mean(double(zooc_150(:,:,y50)),3,'omitnan');
mzmeso = mean(zmeso_150(:,:,y50),3,'omitnan');

%% now look if same
figure
pcolor(glmm); shading flat; colorbar
clim([0 0.25])

figure
pcolor(mzooc); shading flat; colorbar
clim([0 0.25])

figure
pcolor(mzmeso); shading flat; colorbar
clim([0 0.25])

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

load coastlines


%% Interpolate to same grid
% glon=glon';
% glat=glat';
% 
% fZ = griddata(tlat,tlon,mZ_fosi,glat,glon);
% 
% wZ = griddata(CLAT,CLON,mZ_waccm,glat,glon);

%% diffs
diffZooc = (mzooc - glmm);
diffZmeso = (mzmeso - glmm);

zooc_corr = zooc_150 - repmat(diffZooc,1,1,nt);
zmeso_corr = zmeso_150 - repmat(diffZmeso,1,1,nt);

bmzooc = mean(double(zooc_corr(:,:,y50)),3,'omitnan');
bmzmeso = mean(zmeso_corr(:,:,y50),3,'omitnan');

%%
close all

figure(10)
% WACCM zooc
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,mzooc)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('WACCM zooc (molC)');

% WACCM zmeso
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,mzmeso)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
colorbar
set(gcf,'renderer','painters')
title('WACCM zmeso (molC)');

% Diff
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,diffZooc)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.2 0.2]);
colorbar
set(gcf,'renderer','painters')
title('zooc - GLMM');

% PDiff
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,diffZmeso)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.2 0.2]);
colorbar
set(gcf,'renderer','painters')
title('zmeso - GLMM');
%stamp('')
print('-dpng',[ppath 'WACCM_obsGLMM_zooc_zmeso_mean_global_diff_molC.png'])

%% obsGLMM
figure(11)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GLAT,GLON,glmm)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
colorbar
set(gcf,'renderer','painters')
title('GLMM zmeso (molC)');
%stamp('')
print('-dpng',[ppath 'obsGLMM_mean_global_molC.png'])

%% bias corrected
figure(12)
% WACCM zooc
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,bmzooc)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
colorbar
%colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('corrected WACCM zooc (molC)');

% WACCM zmeso
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(CLAT,CLON,bmzmeso)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
colorbar
set(gcf,'renderer','painters')
title('corrected WACCM zmeso (molC)');

% obsglmm
subplot('Position',[0.25 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GLAT,GLON,glmm)
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.1]);
colorbar
set(gcf,'renderer','painters')
title('obsGLMM');

% % PDiff
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(CLAT,CLON,diffZmeso)
% cmocean('balance')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-1 1]);
% colorbar
% set(gcf,'renderer','painters')
% title('zmeso - GLMM');
%stamp('')
print('-dpng',[ppath 'WACCM_obsGLMM_zooc_zmeso_biascorr_global_diff_molC.png'])

%% SAVE
save([wpath 'hist/cesm2_hist_biascorr_zooc_zmeso_diatfraconly_molC_monthly_1850_2014.mat'],...
    'zooc_corr','zmeso_corr','yr','CLON','CLAT','diffZooc','diffZmeso');
