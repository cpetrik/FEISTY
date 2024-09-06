% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss all depths, integrate top 100m
% 1/4 grid from ISIMIP

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% grid - original
ncid = netcdf.open([spath 'full.ocean_static_FishMIP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
geolon = double(geolon);
geolat = double(geolat);

%% map
latlim=[-90 90];
lonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(geolat,geolon,geolat)
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
title('geolat')

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(geolat,geolon,geolon)
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
title('geolon')

figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',mlonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(geolat,geolon,geolon)
h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
title('geolon')









