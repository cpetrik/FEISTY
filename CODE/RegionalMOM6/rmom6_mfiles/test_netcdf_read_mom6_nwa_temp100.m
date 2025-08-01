% Read GFDL MOM6-NWA12
% Temp mean over top 100m

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';

%% one file
ncdisp([fpath 'thetao.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'])

%%
% Dimensions:
% time = 324   (UNLIMITED)
% yh   = 845
% xh   = 775
% nv   = 2

% thetao (#5)
% Size:       775x845x52x324
% Dimensions: xh,yh,z_l,time
% Datatype:   single
% Attributes:
FillValue    = 1.000000020040877e+20;
missing_value = 1.000000020040877e+20;
thetao_units         = 'degC';
thetao_long_name     = 'Sea Water Potential Temperature';
% cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean';
% cell_measures = 'volume: volcello area: areacello';
% time_avg_info = 'average_T1,average_T2,average_DT';
% standard_name = 'sea_water_potential_temperature';

% z_l
% Size:       52x1
% Dimensions: z_l
% Datatype:   double
% Attributes:
z_l_units     = 'meters';
z_l_long_name = 'Depth at cell center';
% axis      = 'Z';
% positive  = 'down';
% edges     = 'z_i';

% time
time_units    = 'days since 1993-01-01 00:00:00';
calendar_type = 'GREGORIAN';
calendar      = 'gregorian';

%%
ncid = netcdf.open([fpath 'thetao.nwa.full.hcast.monthly.raw.r20230520.199301-201912.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:4
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end

for i = 6:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end

%% just top 100
z100 = find(z_l <= 100);

ni = length(xh);
nj = length(yh);
nz = length(z_l);
nt = length(time);

for i = 5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj length(z100) nt]);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

thetao(thetao>1e18) = nan;

%% vis 0m and 100m layers
load([fpath 'nwa_raw_ocean_static_gridspec.mat'],'geolon','geolat','deptho','wet');

[ni,nj]=size(geolon);
geolon = double(geolon);
geolat = double(geolat);
%NWAtl
plotminlat=5; 
plotmaxlat=60;
plotminlon=-100;
plotmaxlon=-30;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure
%subplot(2,2,2)
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,log10(deptho))
cmocean('deep')
load coastlines;                     
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 3]);
hcb = colorbar('h');
title('depth')

t0 = squeeze(thetao(:,:,1,1));
t100 = squeeze(thetao(:,:,15,1));

figure
%subplot(2,2,1)
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,t0)
cmocean('thermal')
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30]);
hcb = colorbar('h');
title('5 m')

figure
%subplot(2,2,3)
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,t100)
cmocean('thermal')
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30]);
hcb = colorbar('h');
title('100 m')

%% grid cell thickness
thk = diff(z_i);
thk100 = single(thk(z100));

thk100m = (nan*ones(ni,nj,15));
for z=1:15
    thk100m(:,:,z) = thk100(z);
end

wet(wet==0) = nan;
wetmat = repmat(wet,1,1,15);
thk100m = thk100m .* wetmat;

thk100m = repmat(thk100m,1,1,1,nt);

thk100m(isnan(thetao(:))) = nan;

%%
test1 = double(squeeze(thk100m(:,:,1,1)));
test15 = double(squeeze(thk100m(:,:,15,1)));

figure
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,test1)
cmocean('deep')
clim([0 10]);
hcb = colorbar('h');
title('thick1')

figure
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,test15)
cmocean('deep')
clim([0 10]);
hcb = colorbar('h');
title('thick100')

%% mean over top 100m
% temp_100 = (nan*ones(ni,nj,nt));
%
% for t=1:nt
%     tp = squeeze(thetao(:,:,:,t));
%     mtp = sum(tp .* thk100m,3,'omitnan') ./ sum(thk100m,3,'omitnan');
%     temp_100(:,:,t) = mtp;
% end

temp_100 = squeeze(sum(thetao .* thk100m,3,'omitnan') ./ sum(thk100m,3,'omitnan'));

temp_100 = double(temp_100);

%%
figure
%subplot(2,2,3)
axesm ('gortho','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat,geolon,squeeze(temp_100(:,:,1)))
cmocean('thermal')
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 30]);
hcb = colorbar('h');
title('t100')

%% Time
yr = 1993 + (time/365);

%%
save([fpath 'temp_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],...
    'FillValue','missing_value','thetao_units','thetao_long_name',...
    'z_l_units','z_l_long_name','time_units','calendar_type','calendar',...
    'time','yr','z_l','temp_100','-v7.3');

