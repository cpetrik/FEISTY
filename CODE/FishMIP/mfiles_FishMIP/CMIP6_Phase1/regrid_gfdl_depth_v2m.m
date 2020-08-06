% Regrid GFDL MOM6 depth (0.5 degree) to 
% Check 1 degree grid from Charlie 
% against CMIP6 grid (1 degree) and MOM6 depth (0.5 degree)

clear all
close all

%% Charlie grid
cpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/';

ncdisp([cpath 'ocean_cobalt_omip_tracers_month_z_1x1deg.static.nc'])

%%
ncid = netcdf.open([cpath 'ocean_cobalt_omip_tracers_month_z_1x1deg.static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

area_units = 'm2';
area_name = 'Ocean Grid-Cell Area';
depth_units = 'm';
depth_name = 'Sea Floor Depth';
sftof_units = '%';
sftof_name = 'Sea Area Fraction';
wet_name = '0 if land, 1 if ocean';

depth_c = double(deptho);
depth_c(depth_c>=1.0e20) = nan;

area_c = areacello;

wet_c = double(wet);
wet_c(wet_c>=1.0e20) = nan;

geolat1 = double(geolat);
geolon1 = double(geolon);
geolat1(geolat1>=1.0e20) = nan;
geolon1(geolon1>=1.0e20) = nan;

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat1,geolon1,depth_c)
colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
title('Charlie depth')

figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat1,geolon1,wet_c)
colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
title('Charlie wet')

%% CMIP6 grid
fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/hist/';

ncdisp([fpath 'gfdl-esm4_r1i1p1f1_historical_expc-bot_onedeg_global_monthly_1850_2014.nc'])

%%
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_expc-bot_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

[lon_c,lat_c] = meshgrid(lon,lat);
[lat_c,lon_c] = meshgrid(lat,lon);

%% MOM6
mpath='/Volumes/FEISTY/GCM_DATA/MOM6/';
ncdisp([mpath 'MOM6_Grid.nc'])
ncid = netcdf.open([mpath 'MOM6_Grid.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
end
netcdf.close(ncid);

deptho(deptho>=1.0e20) = nan;
deptho = double(deptho);
geolon = double(geolon);
geolat = double(geolat);

%% 

%1 CMIP degree grid det_btm
wet1 = expc(:,:,end);
wet1(wet1>=1.0e20) = 0;
wet1(wet1~=0) = 1;


% Check wet
vwet_c = find(wet_c==1);
vwet_c2 = find(wet_c>0);
vwet_1 = find(wet1==1);

%Charlie wet > 0 == CMIP6 variables like det_btm

%% check figs
figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat,geolon,deptho)
colormap('jet')
colorbar
title('MOM6 depth')

figure
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat_c,lon_c,wet1)
colormap('jet')
colorbar
title('CMIP6 wet')


