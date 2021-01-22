% Read GFDL CORE-forced netcdfs

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%% temp
ncdisp([fpath 'ocean.190801-192712.temp.nc'])

% Dimensions:
% xt_ocean       = 360
% yt_ocean       = 200
% time           = 240   (UNLIMITED)
% nv             = 2
% xu_ocean       = 360
% yu_ocean       = 200
% st_ocean       = 50
% st_edges_ocean = 51

% temp
% Size:       360x200x50x240
% Dimensions: xt_ocean,yt_ocean,st_ocean,time
% Datatype:   single
% Attributes:
long_name     = 'Potential temperature';
units         = 'degrees C';
% valid_range   = [-10  500]
% missing_value = -1.000000020040877e+20
% _FillValue    = -1.000000020040877e+20
                       
%% 1908-1927
ncid = netcdf.open([fpath 'ocean.190801-192712.temp.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
z100 = find(st_ocean <= 100);
temp(temp <= -9e19) = NaN;
tp = temp(:,:,z100,:);

% Mean top 100 m
%constant thickness of 10 m in top 200 m
tp1 = squeeze(nanmean(tp,3));
time1 = time;

%% 1928-1947
ncid = netcdf.open([fpath 'ocean.192801-194712.temp.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
z100 = find(st_ocean <= 100);
temp(temp <= -9e19) = NaN;
tp = temp(:,:,z100,:);

% Mean top 100 m
%constant thickness of 10 m in top 200 m
tp2 = squeeze(nanmean(tp,3));
time2 = time;

%% 1948-1967
ncid = netcdf.open([fpath 'ocean.194801-196712.temp.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
z100 = find(st_ocean <= 100);
temp(temp <= -9e19) = NaN;
tp = temp(:,:,z100,:);

% Mean top 100 m
%constant thickness of 10 m in top 200 m
tp3 = squeeze(nanmean(tp,3));
time3 = time;

%% 1968-1987
ncid = netcdf.open([fpath 'ocean.196801-198712.temp.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
z100 = find(st_ocean <= 100);
temp(temp <= -9e19) = NaN;
tp = temp(:,:,z100,:);

% Mean top 100 m
%constant thickness of 10 m in top 200 m
tp4 = squeeze(nanmean(tp,3));
time4 = time;

%% 1988-2007
ncid = netcdf.open([fpath 'ocean.198801-200712.temp.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Top 100 m
z100 = find(st_ocean <= 100);
temp(temp <= -9e19) = NaN;
tp = temp(:,:,z100,:);

% Mean top 100 m
%constant thickness of 10 m in top 200 m
tp5 = squeeze(nanmean(tp,3));
time5 = time;

%% Time
%Put times together

temp1 = cat(3,tp1,tp2);
temp2 = cat(3,temp1,tp3);
temp3 = cat(3,temp2,tp4);
temp_100 = cat(3,temp3,tp5);

time = [time1;time2;time3;time4;time5];

%% Select 1950-2007
yr = (time/365)+1888; %days since 1888-01-01 = first month is Jan 1888
runs = find(yr>=1950);
runs(end)
runs(1)
tp_100 = temp_100(:,:,runs);

%%
save([fpath 'ocean_cobalt_temp100_monthly_1908_2007.mat'],'temp_100','time',...
    'yr','long_name','units','geolat_t','geolon_t');

save([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100','time',...
    'yr','runs','long_name','units','geolat_t','geolon_t');






