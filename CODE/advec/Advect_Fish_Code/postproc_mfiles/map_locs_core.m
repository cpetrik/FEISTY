% See if Hist Locs are in correct place on this grid

clear 

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%%
lpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/';
load([lpath 'core_grid_360x200_id_locs_area_dep.mat']) %,'ids','abbrev');
names = abbrev;

%%
glocs = zeros(ni,nj);
glocs(GRD.ID(ids)) = ones;

%%
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,glocs)
for n=1:length(ids)
    textm(geolat_t(GRD.ID(ids(n))),geolon_t(GRD.ID(ids(n))),abbrev{n});
    hold on;
end
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1]);
hcb = colorbar('h');
print('-dpng',[ppath 'map_CORE_locs.png'])

