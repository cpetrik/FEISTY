%% CORE velocity data

clear
close all

%%
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';

%% oldest file 5/12/2016
ncdisp([vpath 'feb152013_run25_ocean.198801-200712.u.nc'])

% U
% Size:       360x200x50x240
% Dimensions: XU_OCEAN,YU_OCEAN,ST_OCEAN,TIME
% Datatype:   single
% Attributes:
% missing_value = -1.000000020040877e+20
% FillValue    = -1.000000020040877e+20
% long_name     = 'i-current'
% units         = 'm/sec'
% history       = 'From ocean.198801-200712.u'

ncdisp([vpath 'feb152013_run25_ocean.198801-200712.v.nc'])

% V
% Size:       360x200x50x240
% Dimensions: XU_OCEAN,YU_OCEAN,ST_OCEAN,TIME
% Datatype:   single
% Attributes:
% missing_value = -1.000000020040877e+20
% FillValue    = -1.000000020040877e+20
% long_name     = 'j-current'
% units         = 'm/sec'
% history       = 'From ocean.198801-200712.v'

% Velocities at all 50 depth levels
 
%% next file 5/14/2016
ncdisp([vpath 'ocean_cobalt.194801-200712_u_200.nc'])

% U200
% Size:       []
% Dimensions: XU_OCEAN,YU_OCEAN,ST_OCEAN1_20,TIME
% Datatype:   double
% Attributes:
% missing_value = -9.999999999999999e+33
% _FillValue    = -9.999999999999999e+33
% long_name     = 'U[K=1:20]'
% history       = 'From ocean.194801-200712.u'

% ST_OCEAN1_20
% Size:       20x1
% Dimensions: ST_OCEAN1_20
% Datatype:   double
% Attributes:
% long_name     = 'tcell zstar depth'
% units         = 'meters'
% positive      = 'down'
% point_spacing = 'uneven'
% axis          = 'Z'
% standard_name = 'depth'
% bounds        = 'ST_OCEAN1_20_bnds'


ncdisp([vpath 'ocean_cobalt.194801-200712_v_200.nc'])

% V200
% Size:       360x200x20x720
% Dimensions: XU_OCEAN,YU_OCEAN,ST_OCEAN1_20,TIME
% Datatype:   double
% Attributes:
% missing_value = -9.999999999999999e+33
% _FillValue    = -9.999999999999999e+33
% long_name     = 'V[K=1:20]'
% history       = 'From ocean.194801-200712.v'

% Velocities at top 20 depth levels = 0-200 m

%% Next file 5/16/2016

ncdisp([vpath 'feb152013_run25_ocean.198801-200712_u200_v200.nc'])

% Ut_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vt_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% 2D Velocities over 200 m, means or sums?
% ~ -1 to 1 m/s, so means

ncid = netcdf.open([vpath 'feb152013_run25_ocean.198801-200712_u200_v200.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Ut_200(Ut_200<-900)=nan;
Vt_200(Vt_200<-900)=nan;
min(Ut_200(:))
max(Ut_200(:))
min(Vt_200(:))
max(Vt_200(:))


%% Next file 6/17/2016

ncdisp([vpath 'feb152013_run25_ocean.198801-200712_uh200_vh200.nc'])

% Ut_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vt_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Uth_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vth_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% time
% Size:       240x1
% Dimensions: ntime
% Datatype:   double

% geolon_t
% Size:       360x200
% Dimensions: nlon,nlat
% Datatype:   double

% geolat_t
% Size:       360x200
% Dimensions: nlon,nlat
% Datatype:   double

% st_ocean
% Size:       50x1
% Dimensions: ndepth
% Datatype:   double

% 2D Velocities over 200 m, means (Ut) and sums (Uth)?
% Ut and Vt ~ +/- 1, so means
% Uth and Vth ~ +/- 180, so transports (sums)

ncid = netcdf.open([vpath 'feb152013_run25_ocean.198801-200712_uh200_vh200.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

min(Ut_200(:))
max(Ut_200(:))
min(Vt_200(:))
max(Vt_200(:))
min(Uth_200(:))
max(Uth_200(:))
min(Vth_200(:))
max(Vth_200(:))

%% Next file 1/25/2017

ncdisp([vpath 'ocean.194801-200712.u50.nc'])

% U50
% Size:       92x30x5x720
% Dimensions: XU_OCEAN205_296,YU_OCEAN149_178,ST_OCEAN1_5,TIME
% Datatype:   double
% Attributes:
% missing_value = -9.999999999999999e+33
% _FillValue    = -9.999999999999999e+33
% long_name     = 'U[I=205:296,J=149:178,K=1:5]'
% history       = 'From ocean.194801-196712.u'

ncdisp([vpath 'ocean.194801-200712.v50.nc'])

% V50
% Size:       92x30x5x720
% Dimensions: XU_OCEAN205_296,YU_OCEAN149_178,ST_OCEAN1_5,TIME
% Datatype:   double
% Attributes:
% missing_value = -9.999999999999999e+33
% _FillValue    = -9.999999999999999e+33
% long_name     = 'V[I=205:296,J=149:178,K=1:5]'
% history       = 'From ocean.194801-196712.v'

% Velocities at top 5 depth levels = 0-50 m
% ~ +/- 0.5 m/s

ncid = netcdf.open([vpath 'ocean.194801-200712.u50.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

U50(U50<-1e30)=nan;
max(U50(:))
min(U50(:))

%% Next file 1/25/2017

ncdisp([vpath 'ocean.194801-200712_uh50_vh50_v2.nc'])

% Ut_50
% Size:       92x30x720
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vt_50
% Size:       92x30x720
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Uth_50
% Size:       92x30x720
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vth_50
% Size:       92x30x720
% Dimensions: nlon,nlat,ntime
% Datatype:   double

ncid = netcdf.open([vpath 'ocean.194801-200712_uh50_vh50_v2.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

% 2D Velocities over 50 m
% Ut and Vt ~ +/- 0.5, so means
% Uth and Vth ~ +/- 20, so transports (sums)

min(Ut_50(:))
max(Ut_50(:))
min(Vt_50(:))
max(Vt_50(:))
min(Uth_50(:))
max(Uth_50(:))
min(Vth_50(:))
max(Vth_50(:))

%% Last file 2/16/2017

ncdisp([vpath 'feb152013_run25_ocean.198801-200712_uh200_vh200_v2.nc'])

% Ut_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vt_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Uth_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

% Vth_200
% Size:       360x200x240
% Dimensions: nlon,nlat,ntime
% Datatype:   double

ncid = netcdf.open([vpath 'feb152013_run25_ocean.198801-200712_uh200_vh200_v2.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

% 2D Velocities over 50 m
% Ut and Vt ~ +/- 1, so means
% Uth and Vth ~ +/- 200, so transports (sums)

min(Ut_200(:))
min(Vt_200(:))
min(Uth_200(:))
min(Vth_200(:))

max(Ut_200(:))
max(Vt_200(:))
max(Uth_200(:))
max(Vth_200(:))

%%

ncdisp([vpath 'ocean_cobalt_zoo_diags.198801-200712.jhploss_n_Mdz.nc'])
