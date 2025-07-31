% CORE velocity data
% Take mean over top 100 m

clear
close all

%%
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
 
%% 
ncdisp([vpath 'feb152013_run25_ocean.198801-200712.u.nc'])

% U  (#8)          
% Size:       360x200x50x240
% Dimensions: XU_OCEAN,YU_OCEAN,ST_OCEAN,TIME
% Datatype:   single
% Attributes:
missing_value = -1.000000020040877e+20;
FillValue    = -1.000000020040877e+20;
U_long_name     = 'i-current';
U_units         = 'm/sec';
% history       = 'From ocean.198801-200712.u'

% ST_OCEAN     
% Size:       50x1
% Dimensions: ST_OCEAN
% Datatype:   double
% Attributes:
ST_OCEAN_long_name     = 'tcell zstar depth';
ST_OCEAN_units         = 'meters';
% positive      = 'down'
% point_spacing = 'uneven'
% axis          = 'Z'
ST_OCEAN_standard_name = 'depth';
% bounds        = 'ST_OCEAN_bnds'

ncdisp([vpath 'feb152013_run25_ocean.198801-200712.v.nc'])

% V (#8)
% Size:       360x200x50x240
% Dimensions: XU_OCEAN,YU_OCEAN,ST_OCEAN,TIME
% Datatype:   single
% Attributes:
% missing_value = -1.000000020040877e+20
% _FillValue    = -1.000000020040877e+20
V_long_name     = 'j-current';
V_units         = 'm/sec';
% history       = 'From ocean.198801-200712.v'

% TIME
time_units   = 'days since 1888-01-01 00:00:00';
calendar     = 'NOLEAP';

% Velocities at top 50 depth levels 

%% everything except Uvel first 
ncid = netcdf.open([vpath 'feb152013_run25_ocean.198801-200712.u.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:7
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end

for i = 9:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end

%% subset top 100 m
z100 = find(ST_OCEAN <= 100);

[ni,nj] = size(GEOLAT_C);
nk = length(z100);
nt = length(TIME);

%% grid cell thickness
thk = ST_OCEAN_bnds(2,:) - ST_OCEAN_bnds(1,:);
thk100 = single(thk(z100));

thk100m = (nan*ones(ni,nj,nk));
for z=1:nk
    thk100m(:,:,z) = thk100(z);
end

thk100m = repmat(thk100m,1,1,1,nt);

%% U vel =====================================
for i = 8
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nk nt]);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

%% get rid of nans first
U(U<-1e19) = nan;

%% mean over top 100m
u_100 = squeeze(sum(U .* thk100m,3,'omitnan') ./ sum(thk100m,3,'omitnan'));


%% V vel =====================================
ncid = netcdf.open([vpath 'feb152013_run25_ocean.198801-200712.v.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 8
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nk nt]);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

%% get rid of nans first
V(V<-1e19) = nan;

%% mean over top 100m
v_100 = squeeze(sum(V .* thk100m,3,'omitnan') ./ sum(thk100m,3,'omitnan'));

%%
u_100 = double(u_100);
v_100 = double(v_100);

%% Time
yr = 1888 + (TIME/365);


