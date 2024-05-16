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
% units          = 'days since 1959-01-01 00:00:00'
% calendar_type  = 'JULIAN'
% calendar       = 'JULIAN'

% z_l
% Size:       35x1
% long_name      = 'Depth at cell center'
% units          = 'meters'

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

%% Get subset of depth
z100 = find(z_l <= 100);

mthk100 = nan*ones(1,1,length(z100));
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

%% molN/kg --> gWW/m3
nmdz = nmdz * (1/1e-3) * (106/16) * 12.01 * 9;
nlgz = nlgz * (1/1e-3) * (106/16) * 12.01 * 9;

%% Integrate top 100 m
[ni,nj,nk,nt] = size(nmdz);
thkcello = repmat(mthk100,ni,nj,1,nt);
nmdz_100 = squeeze(sum((nmdz.*thkcello),3));
nlgz_100 = squeeze(sum((nlgz.*thkcello),3));

%% zmeso 2nd half
clear nmdz nlgz

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso = netcdf.getVar(ncid,n-1,[0,0,0,run2(1)-1],[1440 720 length(z100) length(run2)]);
end
netcdf.close(ncid);

nmdz(nmdz >= 1.00e+19) = NaN;
nlgz(nlgz >= 1.00e+19) = NaN;

%% molN/kg --> gWW/m3
nmdz = nmdz * (1/1e-3) * (106/16) * 12.01 * 9;
nlgz = nlgz * (1/1e-3) * (106/16) * 12.01 * 9;

%% Integrate top 100 m
[ni,nj,nk,nt] = size(nmdz);
thkcello = repmat(mthk100,ni,nj,1,nt);
nmdz2_100 = squeeze(sum((nmdz.*thkcello),3));
nlgz2_100 = squeeze(sum((nlgz.*thkcello),3));

%% grid
[LAT,LON] = meshgrid(yh,xh);

pcolor(LON,LAT,test2); shading flat

%% Time
yr = 1959 + (time/365);

%%
clear zmeso
nmdz_100(:,:,run2) = nmdz2_100;
nmdz_100 = double(nmdz_100);

nlgz_100(:,:,run2) = nlgz2_100;
nlgz_100 = double(nlgz_100);

mtp = squeeze(nanmean(nanmean(nmdz_100,2),1));
plot(yr,mtp)

%%
clear nmdz nlgz nmdz2_100 nlgz2_100

%%
units_orig = lgz_units;
units_vint = 'gWW m-2';

%%
save([fpath 'ocean_cobalt_mdz_lgz_100m_gWW_month_1961_2010.mat'],...
    'units_vint','time','time_units','yr',...
    'yh','xh','LAT','LON','nmdz_100','nlgz_100','-v7.3');



