% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss all depths, integrate top 100m
% 1/4 grid from ISIMIP

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% Reduce file & data size
% do one time half and then the other
run1 = 1:600;
%run2 = 301:600;


%% zmeso
ncdisp([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'])

%%
% Size:       1440x1080x600
% Dimensions: xh,yh,z_l,time

% Time
% Size:       12x1
% long_name      = 'time'
time_units       = 'days since 1959-01-01 00:00:00';
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'

% xh
% Size:       1440x1
% Dimensions: xh
% Datatype:   double
% Attributes:
xh_long_name      = 'h point nominal longitude';
xh_units          = 'degrees_east';
% cartesian_axis = 'X'

% yh
% Size:       1080x1
% Dimensions: yh
% Datatype:   double
% Attributes:
yh_long_name      = 'h point nominal latitude';
yh_units          = 'degrees_north';
% cartesian_axis = 'Y'

% jhploss_nmdz_100
% Size:       1440x1080x600
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
nmdz_long_name = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m';
nmdz_units     = 'mol m-2 s-1';
missing_value  = 1.000000020040877e+20;
FillValue      = 1.000000020040877e+20;
% cell_methods  = 'area:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'

% jhploss_nlgz_100
% Size:       1440x1080x600
% Dimensions: xh,yh,time
% Datatype:   single
% Attributes:
nlgz_long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m';
nlgz_units         = 'mol m-2 s-1';


%% all but zmeso
ncid = netcdf.open([fpath '19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 3:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Time
yr = 1959 + (time/365);

%% all at once
for n = 1
    varname = netcdf.inqVar(ncid, n-1);
    jhploss_nlgz_100 = netcdf.getVar(ncid,n-1);
    %jhploss_nlgz_100 = netcdf.getVar(ncid,n-1,[0,0,run1(1)-1],[1440 1080 length(run1)]);
end
for n = 2
    varname = netcdf.inqVar(ncid, n-1);
    jhploss_nmdz_100 = netcdf.getVar(ncid,n-1);
    %jhploss_nmdz_100 = netcdf.getVar(ncid,n-1,[0,0,run1(1)-1],[1440 1080 length(run1)]);
end
netcdf.close(ncid);

jhploss_nmdz_100(jhploss_nmdz_100 >= 1.00e+19) = NaN;
jhploss_nlgz_100(jhploss_nlgz_100 >= 1.00e+19) = NaN;

inan = isnan(jhploss_nmdz_100(:));
ocell = ~isnan(squeeze(jhploss_nmdz_100(:,:,1,1)));
lcell = isnan(squeeze(jhploss_nmdz_100(:,:,1,1)));

%% grid
% [LAT,LON] = meshgrid(yh,xh);

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
test = double(squeeze(jhploss_nmdz_100(:,:,6)));
test2 = double(squeeze(jhploss_nlgz_100(:,:,6)));

figure
pcolor(geolon_c,geolat_c,test); shading flat; colorbar; 
colormap('jet')

figure
pcolor(geolon_c,geolat_c,test2); shading flat; colorbar;
colormap('jet')

geolon_c = double(geolon_c);
geolat_c = double(geolat_c);
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',mlonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_c,geolon_c,test)
colormap('jet')
colorbar
%clim([0 0.01])

% figure
% pcolor(LON,LAT,test); shading flat; colorbar;
% colormap('jet')
% 
% figure
% pcolor(LON,LAT,test2); shading flat; colorbar;
% colormap('jet')

%% molN/m2 --> gC/m3
hploss_nmdz_100 = double(jhploss_nmdz_100 * (106/16) * 12.01);
hploss_nlgz_100 = double(jhploss_nlgz_100 * (106/16) * 12.01);

%% viz
test3 = squeeze(hploss_nmdz_100(:,:,6));
test4 = squeeze(hploss_nlgz_100(:,:,6));

figure
pcolor(LON,LAT,3600*test3./(12.01*9)); shading flat; colorbar; %to comp to fishmip molC zmeso-vint
caxis([0.002 0.15])
colormap('jet')

figure
pcolor(LON,LAT,3600*test4./(12.01*9)); shading flat; colorbar;
caxis([0.002 0.15])
colormap('jet')

%% map
clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%%
mtp = squeeze(mean(mean(hploss_nmdz_100,2,'omitnan'),1,'omitnan'));
ltp = squeeze(mean(mean(hploss_nlgz_100,2,'omitnan'),1,'omitnan'));

mzmeso_vint = (mean((hploss_nmdz_100+hploss_nmdz_100),3)); %,'omitnan'));
mz1 = mean(hploss_nmdz_100,3); %,'omitnan');
mz2 = squeeze(hploss_nmdz_100(:,:,6));
lz1 = mean(hploss_nlgz_100,3); %,'omitnan');
lz2 = squeeze(hploss_nlgz_100(:,:,6));

ppath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';

figure(10)
plot(yr,mtp*3600*24,'b'); hold on
plot(yr,ltp*3600*24,'r');

figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,mzmeso_vint*3600*24./12.01)
colormap('jet')
caxis([0 0.01])
colorbar
title('zmeso-vint molC m-2 d-1')
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_HPloss_zmesovint_molC.png'])

figure(11)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mzmeso_vint*3600*24*9))
colormap('jet')
caxis([-2 1])
colorbar
title('zmeso-vint gWW m-2 d-1')
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_HPloss_zmesovint_log10gWW.png'])

figure(2)
pcolor(LON,LAT,log10(mz1*3600*24*9)); shading flat;
cmocean('tempo')
caxis([-3.5 0.5])
colorbar
title('MZ gWW/m2/d')
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_HPloss_MZ_log10gWW.png'])

figure(3)
pcolor(LON,LAT,log10(lz1*3600*24*9)); shading flat;
cmocean('tempo')
caxis([-5 1])
colorbar
title('LZ gWW/m2/d')
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_HPloss_LZ_log10gWW.png'])

figure(4)
pcolor(LON,LAT,log10(mz1*3600*24*1e3)); shading flat;
colormap('jet')
caxis([0 2])
colorbar
title('MZ mgC/m2/d')
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_HPloss_MZ_log10mgC.png'])

figure(5)
pcolor(LON,LAT,log10(lz1*3600*24*1e3)); shading flat;
colormap('jet')
caxis([0 2])
colorbar
title('LZ mgC/m2/d')
print('-dpng',[ppath 'Map_GFDL_15arcmin_jet_HPloss_LZ_log10mgC.png'])


%%
clear jhploss_nmdz_100 jhploss_nlgz_100 

%%
mdz_long_name     = 'medium zooplankton loss to higher preds. integrated in top 100m';
mdz_units         = 'gC m-2 s-1';
lgz_long_name     = 'large zooplankton loss to higher preds.integrated in top 100m';
lgz_units         = 'gC m-2 s-1';

%%
save([fpath 'gfdl-mom6_cobalt2_15arcmin_HPloss_mdz_lgz_100m_month_1961_2010.mat'],...
    'time','time_units','yr',...
    'mdz_long_name','lgz_long_name','mdz_units','lgz_units',...
    'LAT','LON','hploss_nmdz_100','hploss_nlgz_100','-v7.3');



