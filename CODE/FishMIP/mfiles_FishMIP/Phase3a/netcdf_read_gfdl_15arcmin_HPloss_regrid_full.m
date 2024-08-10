% Read GFDL netcdfs
% obsclim
% mdz & lgz HPloss all depths, integrate top 100m
% Regridded by code merging Xiao's with 1 line from Matthias

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
qpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/regular_grid_15arcmin/';

%% zmeso regular grid
ncdisp([fpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc'])

%Global Attributes:
%            CDI              = 'Climate Data Interface version 2.4.2 (https://mpimet.mpg.de/cdi)'
%            Conventions      = 'CF-1.6'
%            filename         = '19610101.ocean_cobalt_fluxes_int.nc'
%            title            = 'OM4p25_LM3TAN_jra_from1959_init234yrs_dynrivernut_O1756_new'
%            associated_files = 'areacello: 19610101.ocean_static.nc'
%            grid_type        = 'regular'
%            grid_tile        = 'N/A'
%            history          = 'Wed Aug 07 16:06:37 2024: cdo remapbil,r1440x720 -selname,jhploss_nlgz_100,jhploss_nmdz_100 /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc
%                               Wed Aug  7 12:20:40 2024: ncatted -O -a coordinates,jhploss_nmdz_100,c,c,geolat geolon /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc
%                               Wed Aug  7 10:56:08 2024: ncatted -O -a coordinates,jhploss_nlgz_100,c,c,geolat geolon /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc
%                               Wed Aug 07 10:13:42 2024: cdo -O -merge /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CPcopy.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/full.ocean_static_FishMIP.nc /Volumes/petrik-lab/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_temp.nc
%                               Fri Nov 17 12:26:44 2023: ncrcat 19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19620101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19630101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19640101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19650101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19660101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19670101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19680101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19690101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19700101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19710101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19720101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19730101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19740101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19750101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19760101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19770101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19780101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19790101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19800101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19810101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19820101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19830101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19840101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19850101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19860101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19870101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19880101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19890101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19900101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19910101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19920101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19930101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19940101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19950101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19960101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19970101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19980101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19990101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20000101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20010101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20020101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20030101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20040101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20050101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20060101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20070101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20080101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20090101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 20100101.ocean_cobalt_fluxes_int_FishMIP_CP.nc 19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP.nc
%                               Fri Nov 17 12:22:25 2023: ncks -v jhploss_nmdz_100,jhploss_nlgz_100 19610101.ocean_cobalt_fluxes_int.nc 19610101.ocean_cobalt_fluxes_int_FishMIP_CP.nc'
%            NCO              = 'netCDF Operators version 5.2.7 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco, Citation = 10.1016/j.envsoft.2008.03.004)'
%            CDO              = 'Climate Data Operators version 2.4.2 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            time = 600   (UNLIMITED)
%            lon  = 1440
%            lat  = 720
% Variables:
%     time            
%            Size:       600x1
%            Dimensions: time
%            Datatype:   double
%            Attributes:
%                        standard_name = 'time'
%                        long_name     = 'time'
time_units         = 'days since 1959-01-01 00:00:00';
%                        calendar      = 'standard'
%                        axis          = 'T'
%     lon             
%            Size:       1440x1
%            Dimensions: lon
%            Datatype:   double
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude'
%                        units         = 'degrees_east'
%                        axis          = 'X'
%     lat             
%            Size:       720x1
%            Dimensions: lat
%            Datatype:   double
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude'
%                        units         = 'degrees_north'
%                        axis          = 'Y'
%     jhploss_nlgz_100
%            Size:       1440x720x600
%            Dimensions: lon,lat,time
%            Datatype:   single
%            Attributes:
nlgz_long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m';
nlgz_units         = 'mol m-2 s-1';
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'
%     jhploss_nmdz_100
%            Size:       1440x720x600
%            Dimensions: lon,lat,time
%            Datatype:   single
%            Attributes:
nmdz_long_name     = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m';
nmdz_units         = 'mol m-2 s-1';
%                        _FillValue    = 1.000000020040877e+20
%                        missing_value = 1.000000020040877e+20
%                        cell_methods  = 'area:mean yh:mean xh:mean time: mean'
%                        cell_measures = 'area: areacello'
%                        time_avg_info = 'average_T1,average_T2,average_DT'

%% 
ncid = netcdf.open([fpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
jhploss_nmdz_100(jhploss_nmdz_100 >= 1.00e+19) = NaN;
jhploss_nlgz_100(jhploss_nlgz_100 >= 1.00e+19) = NaN;

testM = double(squeeze(jhploss_nmdz_100(:,:,60)));
testL = double(squeeze(jhploss_nlgz_100(:,:,160)));

figure
pcolor(testM); shading flat; colorbar; 
%colormap('jet')

figure
pcolor(testL); shading flat; colorbar;

[rLAT,rLON] = meshgrid(lat,lon);

%% Check if same grid as other ISIMIP files =====================
rpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/GFDL_reanalysis/';

%load([rpath 'gridspec_gfdl-mom6-cobalt2_obsclim_15arcmin_orig.mat'])
load([qpath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'])

%% map
latlim=[-90 90];
lonlim=[-280 80];

mlonlim=[-180 180];

load coastlines;

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(rLAT,rLON,testM)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.05]);

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(rLAT,rLON,testL)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 0.05]);

figure
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)
%h=patchm(coastlat,coastlon,'w','FaceColor',[0.75 0.75 0.75]);


%% Need to shift to use with ISIMIP vars, or match up lat-lon
%Xiao: split, flip left and right, and concatenate the matrix before plotting the map

Mtsplit1 = testM(1:720,:);
Mtsplit2 = testM(721:end,:);
Mtflip1 = fliplr(Mtsplit1);
Mtflip2 = fliplr(Mtsplit2);
Mtcomb = [Mtflip2;Mtflip1]; %This looks aligned with ISIMIP depth

Ltsplit1 = testL(1:720,:);
Ltsplit2 = testL(721:end,:);
Ltflip1 = fliplr(Ltsplit1);
Ltflip2 = fliplr(Ltsplit2);
Ltcomb = [Ltflip2;Ltflip1]; %This looks aligned with ISIMIP depth

close all

figure
pcolor(Mtcomb); shading flat

figure
pcolor(Ltcomb); shading flat

figure
pcolor(deptho); shading flat

%% Loop over time to re-orient and save
clear testM testL

[ni,nj,nt] = size(jhploss_nmdz_100);

hploss_nmdz_100 = nan*ones(ni,nj,nt);
hploss_nlgz_100 = nan*ones(ni,nj,nt);

for t=1:nt
    testM = double(squeeze(jhploss_nmdz_100(:,:,t)));
    testL = double(squeeze(jhploss_nlgz_100(:,:,t)));

    Mtsplit1 = testM(1:720,:);
    Mtsplit2 = testM(721:end,:);
    Mtflip1 = fliplr(Mtsplit1);
    Mtflip2 = fliplr(Mtsplit2);
    Mtcomb = [Mtflip2;Mtflip1];

    Ltsplit1 = testL(1:720,:);
    Ltsplit2 = testL(721:end,:);
    Ltflip1 = fliplr(Ltsplit1);
    Ltflip2 = fliplr(Ltsplit2);
    Ltcomb = [Ltflip2;Ltflip1];

    hploss_nmdz_100(:,:,t) = Mtcomb;
    hploss_nlgz_100(:,:,t) = Ltcomb;

    clear Mtsplit1 Mtsplit2 Mtflip1 Mtflip2 Mtcomb
    clear Ltsplit1 Ltsplit2 Ltflip1 Ltflip2 Ltcomb

end

%%
test1 = double(squeeze(hploss_nmdz_100(:,:,500)));
test2 = double(squeeze(hploss_nlgz_100(:,:,360)));

figure
pcolor(test1); shading flat; colorbar; 

figure
pcolor(test2); shading flat; colorbar;

figure
pcolor(deptho); shading flat

%% save
save([fpath '19590101-20101231.ocean_cobalt_fluxes_int_FishMIP_CP_remapped.mat'],...
    'hploss_nmdz_100','hploss_nlgz_100','time','time_units',...
    'nmdz_long_name','nmdz_units','nlgz_long_name','nlgz_units','-v7.3')
