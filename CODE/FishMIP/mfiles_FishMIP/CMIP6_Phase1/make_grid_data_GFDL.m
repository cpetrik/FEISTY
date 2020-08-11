% Make GRD file for FEISTY input from CESM 1 degree model

clear all
close all

Cdir = '/Volumes/FEISTY/Fish-MIP/CMIP6/';

%% Lat & Lon
% ncid = netcdf.open([Cdir 'GFDL/hist/gfdl-esm4_r1i1p1f1_historical_thetao-bot_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
% end
% netcdf.close(ncid);
% thetao(thetao >= 1.00e+20) = NaN;
% 
% %Land mask
% mask = squeeze(thetao(:,:,1));
% lmask = mask;
% lmask(~isnan(lmask)) = 1;
% lmask(isnan(lmask)) = 0;
% 
% %Grid of lat & lon
% [LAT,LON] = meshgrid(lat, lon);

%% Depth, lat, lon
ncid = netcdf.open([Cdir 'GFDL/ocean_cobalt_omip_tracers_month_z_1x1deg.static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
deptho(deptho >= 1.00e+20) = NaN;
geolon(geolon >= 1.00e+20) = NaN;
geolat(geolat >= 1.00e+20) = NaN;
wet(wet >= 1.00e+20) = 0;

%%
geolon_orig = geolon; 
test = geolon; 
id=find(test < -180);
test(id)=test(id)+360;
geolon = double(test);

%% Retain only water cells
ID = find(wet(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = double(geolon(ID));
GRD.LAT = double(geolat(ID));
GRD.Z   = double(deptho(ID));
GRD.lmask = double(wet(ID));
GRD.AREA  = double(areacello(ID));

%% Save needed variables
save([Cdir 'GFDL/gridspec_gfdl.mat'],'deptho','geolat','geolon','wet');
save([Cdir 'GFDL/Data_grid_gfdl.mat'],'GRD');
