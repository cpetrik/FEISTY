% Regrid GFDL MOM6 depth (0.5 degree) to CMIP6 grid (1 degree)

clear all
close all

%% CMIP6 grid
fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/hist/';

ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_historical_expc-bot_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

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

area_units = 'm2';
depth_units = 'm';

deptho(deptho>=1.0e20) = nan;
deptho = double(deptho);
geolon = double(geolon);
geolat = double(geolat);
wet = double(wet);
areacello = double(areacello);

test = geolon; 
id=find(test < -180);
test(id)=test(id)+360;
geolon_shift = double(test);

%% interp
depth_c = griddata(geolon_shift,geolat,deptho,lon_c,lat_c);
area_c = griddata(geolon_shift,geolat,areacello,lon_c,lat_c);

%1 degree grid
wet1 = expc(:,:,end);
wet1(wet1>=1.0e20) = 0;
wet1(wet1~=0) = 1;
wet_c = wet1;

%% save as netcdf
%NaN fills
depth_c(isnan(depth_c))=1e20;
wet_c(isnan(wet_c))=1e20;
area_c(isnan(area_c))=1e20;

[ni,nj]=size(lon_c);

%
file_tsb = '/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/gfdl-esm4_grid.nc4';
ncidSB = netcdf.create(file_tsb,'NETCDF4');

lon_dim = netcdf.defDim(ncidSB,'longitude',ni);
lat_dim = netcdf.defDim(ncidSB,'latitude',nj);

vidlat = netcdf.defVar(ncidSB,'LAT','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');

vidlon = netcdf.defVar(ncidSB,'LON','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );

viddepSB = netcdf.defVar(ncidSB,'depth','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,viddepSB,'long_name','Sea Floor Depth');
netcdf.putAtt(ncidSB,viddepSB,'units','m2' );
netcdf.defVarFill(ncidSB,viddepSB,false,1.0e20);

vidwetSB = netcdf.defVar(ncidSB,'wet','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidwetSB,'long_name','0 if land, 1 if ocean at tracer points');
netcdf.putAtt(ncidSB,vidwetSB,'units','' );
netcdf.defVarFill(ncidSB,vidwetSB,false,1.0e20);

vidareaSB = netcdf.defVar(ncidSB,'area','double',[lon_dim,lat_dim]);
netcdf.putAtt(ncidSB,vidareaSB,'long_name','Ocean Grid-Cell Area');
netcdf.putAtt(ncidSB,vidareaSB,'units','m2' );
netcdf.defVarFill(ncidSB,vidareaSB,false,1.0e20);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,lat_c);
netcdf.putVar(ncidSB,vidlon,lon_c);
netcdf.putVar(ncidSB,viddepSB,depth_c);
netcdf.putVar(ncidSB,vidwetSB,wet_c);
netcdf.putVar(ncidSB,vidareaSB,area_c);

netcdf.close(ncidSB);
