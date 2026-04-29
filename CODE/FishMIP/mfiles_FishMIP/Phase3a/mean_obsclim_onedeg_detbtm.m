% GFDL 1deg obsclim
% Detritus sinking flux btm
% Calc climatol mean 1997-2008

clear 
close all

%%
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/FishMIP/Phase3a/OneDeg/';

load([fpath 'gfdl-mom6-cobalt2_obsclim_expc-bot_onedeg_global_monthly_1961_2010.mat'])

%% subset yrs
yid = find(yr>=1997 & yr<2009);

det = det_btm(:,:,yid);

%% monthly climatol
[ni,nj] = size(LON);

det_mo = nan*ones(ni,nj,12);
for m = 1:12
    mo = m:12:length(yid);
    det_mo(:,:,m) = mean(det(:,:,mo),3,'omitnan');
end

%% annual mean
det_yr = mean(det,3,'omitnan');

%% convert units
% molC/m2/s to mgC/m2/d
det_mo = det_mo * 12.01 * 1e3 * 60*60*24;
det_yr = det_yr * 12.01 * 1e3 * 60*60*24;

%% save
units_new = 'mgC/m2/d';
save([fpath 'gfdl-mom6-cobalt2_obsclim_expc-bot_onedeg_global_means_19997_2008.mat'],...
    'units_new','LAT','LON','det_mo','det_yr');

%% maps
load([fpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([fpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'], 'GRD');

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%%
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_yr))
cmocean('dense')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([-0.5 2.5]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('log10 mean seafloor detritus flux (mgC m^-^2 d^-^1) 1997-2008')
stamp('')
print('-dpng',[ppath 'gfdl-mom6-cobalt2_obsclim_expc-bot_onedeg_global_log10mean_1997_2008.png'])

%%
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(det_yr))
cmocean('dense')
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 250]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
title('mean seafloor detritus flux (mgC m^-^2 d^-^1) 1997-2008')
stamp('')
print('-dpng',[ppath 'gfdl-mom6-cobalt2_obsclim_expc-bot_onedeg_global_mean_1997_2008.png'])