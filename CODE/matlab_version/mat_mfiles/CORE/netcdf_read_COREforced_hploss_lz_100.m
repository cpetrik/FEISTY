% Read GFDL CORE-forced netcdfs

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%%
ncdisp([fpath 'ocean_cobalt_zoo_diags.188801-190712.jhploss_n_Lgz.nc'])

% Dimensions:
% xu_ocean       = 360
% yu_ocean       = 200
% time           = 240   (UNLIMITED)
% nv             = 2
% xt_ocean       = 360
% yt_ocean       = 200
% st_ocean       = 50 - depth middle
% st_edges_ocean = 51 - depth edges

% jhploss_n_Lgz
% Size:       360x200x50x240
% Dimensions: xt_ocean,yt_ocean,st_ocean,time
% Datatype:   single
% Attributes:
long_name     = 'Large zooplankton loss to higher predators layer integral';
units         = 'mol N m-2 s-1';
% missing_value = -10000000000
% _FillValue    = -10000000000

%% 1888-1907
ncid = netcdf.open([fpath 'ocean_cobalt_zoo_diags.188801-190712.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

%% Top 100 m
z100 = find(st_ocean <= 100);
jhploss_n_Lgz(jhploss_n_Lgz <= -9e9) = NaN;
lz = jhploss_n_Lgz(:,:,z100,:);
dz = diff(st_edges_ocean);
thkcello = dz(z100);

% Integrate top 100 m
%each layer already integrated, do not mult by 10
lz1 = squeeze(nansum((lz),3));
time1 = time;

%% 1908-1927
ncid = netcdf.open([fpath 'ocean_cobalt_zoo_diags.190801-192712.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
jhploss_n_Lgz(jhploss_n_Lgz <= -9e9) = NaN;
lz = jhploss_n_Lgz(:,:,z100,:);

% Integrate top 100 m
%each layer already integrated, do not mult by 10
lz2 = squeeze(nansum((lz),3));
time2 = time;

%% 1928-1947
ncid = netcdf.open([fpath 'ocean_cobalt_zoo_diags.192801-194712.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
jhploss_n_Lgz(jhploss_n_Lgz <= -9e9) = NaN;
lz = jhploss_n_Lgz(:,:,z100,:);

% Integrate top 100 m
%each layer already integrated, do not mult by 10
lz3 = squeeze(nansum((lz),3));
time3 = time;

%% 1948-1967
ncid = netcdf.open([fpath 'ocean_cobalt_zoo_diags.194801-196712.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
jhploss_n_Lgz(jhploss_n_Lgz <= -9e9) = NaN;
lz = jhploss_n_Lgz(:,:,z100,:);

% Integrate top 100 m
%each layer already integrated, do not mult by 10
lz4 = squeeze(nansum((lz),3));
time4 = time;

%% 1968-1987
ncid = netcdf.open([fpath 'ocean_cobalt_zoo_diags.196801-198712.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
jhploss_n_Lgz(jhploss_n_Lgz <= -9e9) = NaN;
lz = jhploss_n_Lgz(:,:,z100,:);

% Integrate top 100 m
%each layer already integrated, do not mult by 10
lz5 = squeeze(nansum((lz),3));
time5 = time;

%% 1988-2007
ncid = netcdf.open([fpath 'ocean_cobalt_zoo_diags.198801-200712.jhploss_n_Lgz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
jhploss_n_Lgz(jhploss_n_Lgz <= -9e9) = NaN;
lz = jhploss_n_Lgz(:,:,z100,:);

% Integrate top 100 m
%each layer already integrated, do not mult by 10
lz6 = squeeze(nansum((lz),3));
time6 = time;

%% Time
%Put times together

lgz1 = cat(3,lz1,lz2);
lgz2 = cat(3,lgz1,lz3);
lgz3 = cat(3,lgz2,lz4);
lgz4 = cat(3,lgz3,lz5);
hploss_lgz = cat(3,lgz4,lz6);

time = [time1;time2;time3;time4;time5;time6];

%% Select 1950-2007
yr = (time/365)+1888; %days since 1888-01-01 = first month is Jan 1888
runs = find(yr>=1950);
runs(end)
runs(1)
hploss_lz_100 = hploss_lgz(:,:,runs);

%%
save([fpath 'ocean_cobalt_hploss_lz100_monthly_1888_2007.mat'],'hploss_lgz','time',...
    'yr','long_name','units','geolat_t','geolon_t');

save([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100','time',...
    'yr','runs','long_name','units','geolat_t','geolon_t');






