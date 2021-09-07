% Read GFDL CORE-forced netcdfs

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%%
ncdisp([fpath 'ocean_cobalt_tracers.190801-192712.nmdz.nc'])

% Dimensions:
% xu_ocean       = 360
% yu_ocean       = 200
% time           = 240   (UNLIMITED)
% nv             = 2
% xt_ocean       = 360
% yt_ocean       = 200
% st_ocean       = 50 - depth middle
% st_edges_ocean = 51 - depth edges

% nmdz
% Size:       360x200x50x240
% Dimensions: xt_ocean,yt_ocean,st_ocean,time
% Datatype:   single
% Attributes:
long_name     = 'Medium-sized zooplankton Nitrogen';
units         = 'mol/kg';
% missing_value = -10000000000
% _FillValue    = -10000000000

%% 1908-1927
ncid = netcdf.open([fpath 'ocean_cobalt_tracers.190801-192712.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

%% Top 100 m
z100 = find(st_ocean <= 100);
nmdz(nmdz <= -9e9) = NaN;
mz = nmdz(:,:,z100,:);
dz = diff(st_edges_ocean);
thkcello = dz(z100);

% Integrate top 100 m
%uses constant thickness of 10 m in top 200 m
%kg to m3: 10^3 kg in 1 m^3
%m3 to m2: mult by bin depth thickness of 10
mz1 = squeeze(nansum((1e3*10*mz),3));
time1 = time;

%% 1928-1947
ncid = netcdf.open([fpath 'ocean_cobalt_tracers.192801-194712.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
nmdz(nmdz <= -9e9) = NaN;
mz = nmdz(:,:,z100,:);

% Integrate top 100 m
%uses constant thickness of 10 m in top 200 m
%kg to m3: 10^3 kg in 1 m^3
%m3 to m2: mult by bin depth thickness of 10
mz2 = squeeze(nansum((1e3*10*mz),3));
time2 = time;

%% 1948-1967
ncid = netcdf.open([fpath 'ocean_cobalt_tracers.194801-196712.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
nmdz(nmdz <= -9e9) = NaN;
mz = nmdz(:,:,z100,:);

% Integrate top 100 m
%uses constant thickness of 10 m in top 200 m
mz3 = squeeze(nansum((1e3*10*mz),3));
time3 = time;

%% 1968-1987
ncid = netcdf.open([fpath 'ocean_cobalt_tracers.196801-198712.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
nmdz(nmdz <= -9e9) = NaN;
mz = nmdz(:,:,z100,:);

% Integrate top 100 m
%uses constant thickness of 10 m in top 200 m
mz4 = squeeze(nansum((1e3*10*mz),3));
time4 = time;

%% 1988-2007
ncid = netcdf.open([fpath 'ocean_cobalt_tracers.198801-200712.nmdz.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars-1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
nmdz(nmdz <= -9e9) = NaN;
mz = nmdz(:,:,z100,:);

% Integrate top 100 m
%uses constant thickness of 10 m in top 200 m
mz5 = squeeze(nansum((1e3*10*mz),3));
time5 = time;

%%
units_orig = units;
units_vint = 'mol m-2';

%% Time
%Put times together

mdz1 = cat(3,mz1,mz2);
mdz2 = cat(3,mdz1,mz3);
mdz3 = cat(3,mdz2,mz4);
mdz = cat(3,mdz3,mz5);

time = [time1;time2;time3;time4;time5];

%% Select 1950-2007
yr = (time/365)+1888; %days since 1888-01-01 = first month is Jan 1888
runs = find(yr>=1950);
runs(end)
runs(1)
mz_100 = mdz(:,:,runs);

%%
save([fpath 'ocean_cobalt_mz100_monthly_1908_2007.mat'],'mdz','time',...
    'yr','long_name','units_orig','units_vint','geolat_t','geolon_t');

save([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100','time',...
    'yr','runs','long_name','units_orig','units_vint','geolat_t','geolon_t');






