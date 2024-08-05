% Read GFDL netcdfs
% obsclim
% mdz & lgz biomass all depths
% Regridded by code merging Xiao's with 1 line from Matthias

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% zmeso regular grid
ncdisp([fpath '19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc'])

%Global Attributes:

%% 
ncid = netcdf.open([fpath '19710101.ocean_cobalt_tracers_month_z_FishMIP_CP_remapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
nmdz(nmdz >= 1.00e+19) = NaN;

inan = isnan(nmdz(:));
ocell = ~isnan(squeeze(nmdz(:,:,1,1)));
lcell = isnan(squeeze(nmdz(:,:,1,1)));

%% molN/m3 --> gC/m2
nmdz_100 = double(nmdz * (106/16) * 12.01 * 100);

%%
test = (squeeze(nmdz_100(:,:,1,6)));

figure
pcolor(test); shading flat; colorbar; 
colormap('jet')

[rLAT,rLON] = meshgrid(lat,lon);

%% Check if same grid as other ISIMIP files =====================
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
load([qpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])
load([qpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])

%% Same orientation?
figure(1)
pcolor(rLAT); shading flat;
colorbar
title('rLAT')

figure(2)
pcolor(LAT); shading flat;
colorbar
title('LAT')

figure(3)
pcolor(rLON); shading flat;
colorbar
title('LON')

figure(4)
pcolor(LON); shading flat;
colorbar
title('LON')

%flipped LR and shifted in matrix
figure(5)
pcolor(fliplr(test)); shading flat; 
colorbar
clim([0 0.01])
title('HP')

figure(6)
pcolor(deptho); shading flat
colorbar
title('depth')

%% Depth has a few more ocean cells
oid = find(~isnan(squeeze(nmdz(:,:,1,1))));

% oid       670568x1
% GRD.ID    670589x1
gdiff=setdiff(GRD.ID,oid); %off because shifted

%% map
latlim=[-90 90];
lonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(rLAT,rLON,test)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.05]);

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);



%% Need to shift to use with ISIMIP vars, or match up lat-lon



