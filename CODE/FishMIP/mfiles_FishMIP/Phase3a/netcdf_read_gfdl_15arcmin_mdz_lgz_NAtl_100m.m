% Read GFDL netcdfs
% obsclim
% mdz & lgz all depths, integrate top 100m
% 1/4 grid from ISIMIP

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
cpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';

%% Reduce file & data size
% do one time half and then the other
run1 = 1:300;
run2 = 301:600;

% just N Atl
latmin = 540;
latmax = 1080;
lonmin = 720;
lonmax = 1440;

latcount = latmax - latmin;
loncount = lonmax - lonmin;

%% zmeso
ncdisp([fpath '19610101-20101231.ocean_cobalt_tracers_month_z_FishMIP_CP.nc'])

%%
% Size:       1440x1080x35x12
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

% z_l
% Size:       35x1
zl_long_name      = 'Depth at cell center';
zl_units          = 'meters';

% nmdz
% Size:       1440x1080x35x12
% Dimensions: xh,yh,z_l,time
mdz_long_name     = 'Medium-sized zooplankton Nitrogen';
mdz_units         = 'mol/kg';
missing_value = 1.000000020040877e+20;
FillValue    = 1.000000020040877e+20;
% cell_methods  = 'area:mean z_l:mean yh:mean xh:mean time: mean'
% cell_measures = 'area: areacello'
% time_avg_info = 'average_T1,average_T2,average_DT'

% nlgz
% Size:       1440x1080x35x12
% Dimensions: xh,yh,z_l,time
lgz_long_name     = 'large Zooplankton Nitrogen';
lgz_units         = 'mol/kg';

%% Cell thickness = thkcello
tcid = netcdf.open([cpath 'gfdl-mom6-cobalt2_obsclim_thkcello_15arcmin_global_fixed.nc'],'NC_NOWRITE');
[tdims,tvars,tgatts,unlimdimid] = netcdf.inq(tcid);

% just last var = thkcello
for t = tvars
    varname = netcdf.inqVar(tcid, t-1);
    eval([ varname ' = netcdf.getVar(tcid,t-1);']);
    %thkcello = netcdf.getVar(tcid,t-1, [0,0,0,runs(1)-1],[360 180 length(z100) length(runs)]);
end
netcdf.close(tcid);
thkcello(thkcello >= 1.00e+20) = NaN;
thkcello = double(thkcello);

thk = reshape(thkcello,1440*720,35);
mthk = mean(thk,1,'omitnan');
mthk = mthk';

%%
clear thkcello thk

%% all but zmeso
ncid = netcdf.open([fpath '19610101-20101231.ocean_cobalt_tracers_month_z_FishMIP_CP.nc'],'NC_NOWRITE');
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
aLAT = LAT(lonmin+1:lonmax,latmin+1:latmax);
aLON = LON(lonmin+1:lonmax,latmin+1:latmax);

%% Get subset of depth
z100 = find(z_l <= 100);

mthk100 = ones(1,1,length(z100));
mthk100(:,:,1:end) = mthk(z100);

%% zmeso 1st half
for n = 1
    varname = netcdf.inqVar(ncid, n-1);
    nlgz = netcdf.getVar(ncid,n-1,[lonmin,latmin,0,run1(1)-1],[loncount latcount length(z100) length(run1)]);
end
for n = 2
    varname = netcdf.inqVar(ncid, n-1);
    nmdz = netcdf.getVar(ncid,n-1,[lonmin,latmin,0,run1(1)-1],[loncount latcount length(z100) length(run1)]);
end
nmdz(nmdz >= 1.00e+19) = NaN;
nlgz(nlgz >= 1.00e+19) = NaN;

inan = isnan(nmdz(:));
ocell = ~isnan(squeeze(nmdz(:,:,1,1)));
lcell = isnan(squeeze(nmdz(:,:,1,1)));

%%
test = double(squeeze(nmdz(:,:,1,6)));
test2 = double(squeeze(nlgz(:,:,1,6)));

figure
pcolor(aLON,aLAT,test); shading flat; colorbar; 
colormap('jet')

figure
pcolor(aLON,aLAT,test2); shading flat; colorbar;
colormap('jet')

%% molN/kg --> gWW/m3
nmdz = double(nmdz * (1/1e-3) * (106/16) * 12.01 * 9);
nlgz = double(nlgz * (1/1e-3) * (106/16) * 12.01 * 9);

%% Integrate top 100 m
[ni,nj,nk,nt] = size(nmdz);

thkcello = repmat(mthk100,ni,nj,1,nt);
thkcello(inan) = NaN;

nmdz_m = nmdz.*thkcello;
nlgz_m = nlgz.*thkcello;

nmdz_100 = squeeze(sum(nmdz_m,3,'omitnan'));
nlgz_100 = squeeze(sum(nlgz_m,3,'omitnan'));

znan = repmat(lcell,1,1,nt);
nmdz_100(znan) = NaN;
nlgz_100(znan) = NaN;

%% viz
test3 = squeeze(nmdz_100(:,:,6));
test4 = squeeze(nlgz_100(:,:,6));

figure
pcolor(aLON,aLAT,test3./(12.01*9)); shading flat; colorbar; %div by 12.01 to comp to fishmip molC zmeso-vint
caxis([0.002 0.15])
colormap('jet')

figure
pcolor(aLON,aLAT,test4./(12.01*9)); shading flat; colorbar;
caxis([0.002 0.15])
colormap('jet')

%% zmeso 2nd half -------------------------------------------------------
clear nmdz nlgz thkcello nmdz_m nlgz_m znan inan lcell

for n = 1
    varname = netcdf.inqVar(ncid, n-1);
    nlgz = netcdf.getVar(ncid,n-1,[lonmin,latmin,0,run2(1)-1],[loncount latcount length(z100) length(run2)]);
end
for n = 2
    varname = netcdf.inqVar(ncid, n-1);
    nmdz = netcdf.getVar(ncid,n-1,[lonmin,latmin,0,run2(1)-1],[loncount latcount length(z100) length(run2)]);
end

netcdf.close(ncid);

nmdz(nmdz >= 1.00e+19) = NaN;
nlgz(nlgz >= 1.00e+19) = NaN;

inan = isnan(nmdz(:));
lcell = isnan(squeeze(nmdz(:,:,1,1)));

%% molN/kg --> gWW/m3
nmdz = double(nmdz * (1/1e-3) * (106/16) * 12.01 * 9);
nlgz = double(nlgz * (1/1e-3) * (106/16) * 12.01 * 9);

%% Integrate top 100 m
[ni,nj,nk,nt] = size(nmdz);

thkcello = repmat(mthk100,ni,nj,1,nt);
thkcello(inan) = NaN;

nmdz2_m = nmdz.*thkcello;
nlgz2_m = nlgz.*thkcello;

nmdz2_100 = squeeze(sum(nmdz2_m,3,'omitnan'));
nlgz2_100 = squeeze(sum(nlgz2_m,3,'omitnan'));

znan = repmat(lcell,1,1,nt);
nmdz2_100(znan) = NaN;
nlgz2_100(znan) = NaN;

%% Put together  -------------------------------------------------------
nmdz_100(:,:,run2) = nmdz2_100;
nmdz_100 = double(nmdz_100);

nlgz_100(:,:,run2) = nlgz2_100;
nlgz_100 = double(nlgz_100);

%% map
clatlim=[-90 90];
clonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

%%
mtp = squeeze(mean(mean(nmdz_100,2,'omitnan'),1,'omitnan'));

mzmeso_vint = (mean((nmdz_100+nlgz_100),3)); %,'omitnan'));
test1 = mean(nmdz2_100,3); %,'omitnan');
test2 = squeeze(nmdz2_100(:,:,6));
test3 = squeeze(nmdz2_m(:,:,1,6));
test4 = squeeze(nmdz(:,:,1,6));
test5 = squeeze(thkcello(:,:,1,6));

ppath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/FishMIP/Phase3a/';

figure(10)
plot(yr,mtp)

figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(aLAT,aLON,mzmeso_vint./12.01./9)
colormap('jet')
caxis([0.02 0.15])
colorbar
title('zmeso-vint molC')
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_NAtl_GFDL_15arcmin_jet_zmesovint_molC.png'])

figure(11)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(aLAT,aLON,log10(mzmeso_vint))
colormap('jet')
caxis([-1 1.5])
colorbar
title('zmeso-vint molC')
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[ppath 'Map_NAtl_GFDL_15arcmin_jet_zmesovint_log10gWW.png'])

figure(2)
pcolor(aLON,aLAT,mzmeso_vint); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('zmeso-vint molC')

figure(3)
pcolor(aLON,aLAT,test1); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('test1')

figure(4)
pcolor(aLON,aLAT,test2); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('test2')

figure(5)
pcolor(aLON,aLAT,test3); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('test3')

figure(6)
pcolor(aLON,aLAT,test4); shading flat;
colormap('jet')
caxis([0 10])
colorbar
title('test4')

figure(7)
pcolor(aLON,aLAT,test5); shading flat;
%colormap('jet')
%caxis([0 10])
colorbar
title('thk')

%%
clear nmdz nlgz nmdz2_100 nlgz2_100

%%
mdz_long_name     = 'medium zooplankton biomass integrated in top 100m';
mdz_units         = 'gWW m-2';
lgz_long_name     = 'large zooplankton biomass integrated in top 100m';
lgz_units         = 'gWW m-2';

%%
save([fpath 'mom6_cobalt2_mdz_lgz_100m_gWW_month_1961_2010.mat'],...
    'time','time_units','yr',...
    'mdz_long_name','lgz_long_name','mdz_units','lgz_units',...
    'aLAT','aLON','nmdz_100','nlgz_100','-v7.3');



