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
%% grid cell thickness
thk = diff(z_i);
thk100 = single(thk(z100));

thk100m = (nan*ones(ni,nj,15));
for z=1:15
    thk100m(:,:,z) = thk100(z);
end

%% mean over top 100m
temp_100 = (nan*ones(ni,nj,nt));
for t=1:nt
    tp = squeeze(thetao(:,:,:,t));
    mtp = sum(tp .* thk100m,3,'omitnan') ./ sum(thk100m,3,'omitnan');
    temp_100(:,:,t) = mtp;
end

%%
temp_100(temp_100>1e19) = nan;
temp_100 = double(temp_100);

%% Time
yr = 1993 + (time/365);

%%
save([fpath 'temp_100.nwa.full.hcast.monthly.raw.r20230520.199301-201912.mat'],...
    'FillValue','missing_value','thetao_units','thetao_long_name',...
    'z_l_units','z_l_long_name','time_units','calendar_type','calendar',...
    'time','yr','z_l','temp_100','-v7.3');

