% Read GFDL CORE-forced netcdfs
% bottom POC flux

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%%
ncdisp([fpath 'ocean_cobalt_misc.188801-190712.fndet_btm.nc'])

% Dimensions:
% xu_ocean = 360
% yu_ocean = 200
% time     = 240   (UNLIMITED)
% nv       = 2
% xt_ocean = 360
% yt_ocean = 200
% time:units = 'days since 1888-01-01 00:00:00'

% fndet_btm
% Size:       360x200x240
% Dimensions: xt_ocean,yt_ocean,time
% Datatype:   single
% Attributes:
long_name     = 'ndet sinking flux to bottom';
units         = 'mol m-2 s-1';
% missing_value = -10000000000
% _FillValue    = -10000000000

%% 1888-1907
ncid = netcdf.open([fpath 'ocean_cobalt_misc.188801-190712.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% NaNs on land cells
fndet_btm(fndet_btm <= -9e9) = NaN;
det1 = fndet_btm;
time1 = time;

%% 1908-1927
ncid = netcdf.open([fpath 'ocean_cobalt_misc.190801-192712.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% NaNs on land cells
fndet_btm(fndet_btm <= -9e9) = NaN;
det2 = fndet_btm;
time2 = time;

%% 1928-1947
ncid = netcdf.open([fpath 'ocean_cobalt_misc.192801-194712.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% NaNs on land cells
fndet_btm(fndet_btm <= -9e9) = NaN;
det3 = fndet_btm;
time3 = time;

%% 1948-1967
ncid = netcdf.open([fpath 'ocean_cobalt_misc.194801-196712.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% NaNs on land cells
fndet_btm(fndet_btm <= -9e9) = NaN;
det4 = fndet_btm;
time4 = time;

%% 1968-1987
ncid = netcdf.open([fpath 'ocean_cobalt_misc.196801-198712.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% NaNs on land cells
fndet_btm(fndet_btm <= -9e9) = NaN;
det5 = fndet_btm;
time5 = time;

%% 1988-2007
ncid = netcdf.open([fpath 'ocean_cobalt_misc.198801-200712.fndet_btm.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
end
netcdf.close(ncid);

% NaNs on land cells
fndet_btm(fndet_btm <= -9e9) = NaN;
det6 = fndet_btm;
time6 = time;

%% Time
%Put times together

expc1 = cat(3,det1,det2);
expc2 = cat(3,expc1,det3);
expc3 = cat(3,expc2,det4);
expc4 = cat(3,expc3,det5);
expc = cat(3,expc4,det6);

time = [time1;time2;time3;time4;time5;time6];

%% Select 1950-2007
yr = (time/365)+1888; %days since 1888-01-01 = first month is Jan 1888

runs = find(yr>=1950);

runs(end)
runs(1)
det_btm = expc(:,:,runs);

%%
save([fpath 'ocean_cobalt_fndet_btm_monthly_1888_2007.mat'],'expc','time',...
    'yr','long_name','units','geolat_t','geolon_t');

save([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat'],'det_btm','time',...
    'yr','runs','long_name','units','geolat_t','geolon_t');





