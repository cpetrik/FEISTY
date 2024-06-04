% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss all depths, integrate top 100m
% 1/4 grid from ISIMIP

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

%% Reduce file & data size
% do one time half and then the other
run1 = 1:600;
%run2 = 301:600;


%% zmeso
ncdisp([fpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP.nc'])

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
ncid = netcdf.open([fpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 3:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Time
yr = 1959 + (time/365);

%% grid
[LAT,LON] = meshgrid(yh,xh);

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

%%
test = double(squeeze(jhploss_nmdz_100(:,:,6)));
test2 = double(squeeze(jhploss_nlgz_100(:,:,6)));

figure
pcolor(LON,LAT,test); shading flat; colorbar; 
colormap('jet')

figure
pcolor(LON,LAT,test2); shading flat; colorbar;
colormap('jet')

%% molN/kg --> gC/m3
hploss_nmdz_100 = double(jhploss_nmdz_100 * (1/1e-3) * (106/16) * 12.01);
hploss_nlgz_100 = double(jhploss_nlgz_100 * (1/1e-3) * (106/16) * 12.01);

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

mzmeso_vint = (mean((hploss_nmdz_100+hploss_nmdz_100),3)); %,'omitnan'));
test1 = mean(hploss_nmdz_100,3); %,'omitnan');
test2 = squeeze(hploss_nmdz_100(:,:,6));

ppath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';

figure(10)
plot(yr,mtp*3600)

figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,mzmeso_vint*3600./12.01)
colormap('jet')
caxis([0.02 0.15])
colorbar
title('zmeso-vint molC d-1')
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_NAtl_GFDL_15arcmin_jet_HPloss_zmesovint_molC.png'])

figure(11)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mzmeso_vint*3600))
colormap('jet')
caxis([-1 1.5])
colorbar
title('zmeso-vint log_1_0 gC d-1')
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_NAtl_GFDL_15arcmin_jet_HPloss_zmesovint_log10gC.png'])

figure(2)
pcolor(LON,LAT,mzmeso_vint*3600); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('zmeso-vint molC')

figure(3)
pcolor(LON,LAT,test1*3600); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('test1')

figure(4)
pcolor(LON,LAT,test2*3600); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('test2')


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



