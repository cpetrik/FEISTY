% Read GFDL CORE-forced netcdfs
% bottom temp
% need to find depth cell that is "bottom" for each lat,lon

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';
load([fpath 'ocean_cobalt_grid.mat'])

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
temp(temp <= -9e19) = NaN;

%% Seafloor
[ni,nj,nz,nt] = size(temp);

%Works until z=44, then arrays too big
% tb1 = NaN*ones([ni,nj,nt]);
% for z = 1:50
%     id = find(kmt(:)==z);
%     [m,n] = ind2sub([ni,nj],id);
%     tb1(m,n,:) = temp(m,n,z,:);
% end

tzall = reshape(temp,ni*nj,nz,nt);
tb = NaN*ones([ni*nj,nt]);
for z = 1:50
    id = (kmt(:)==z);
    tb(id,:) = tzall(id,z,:);
end

tb1 = reshape(tb,ni,nj,nt);
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
temp(temp <= -9e19) = NaN;

% Seafloor
[ni,nj,nz,nt] = size(temp);
tzall = reshape(temp,ni*nj,nz,nt);
tb = NaN*ones([ni*nj,nt]);
for z = 1:50
    id = (kmt(:)==z);
    tb(id,:) = tzall(id,z,:);
end

tb2 = reshape(tb,ni,nj,nt);
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
temp(temp <= -9e19) = NaN;

% Seafloor
[ni,nj,nz,nt] = size(temp);
tzall = reshape(temp,ni*nj,nz,nt);
tb = NaN*ones([ni*nj,nt]);
for z = 1:50
    id = (kmt(:)==z);
    tb(id,:) = tzall(id,z,:);
end

tb3 = reshape(tb,ni,nj,nt);
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
temp(temp <= -9e19) = NaN;

% Seafloor
[ni,nj,nz,nt] = size(temp);
tzall = reshape(temp,ni*nj,nz,nt);
tb = NaN*ones([ni*nj,nt]);
for z = 1:50
    id = (kmt(:)==z);
    tb(id,:) = tzall(id,z,:);
end

tb4 = reshape(tb,ni,nj,nt);
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
temp(temp <= -9e19) = NaN;

% Seafloor
[ni,nj,nz,nt] = size(temp);
tzall = reshape(temp,ni*nj,nz,nt);
tb = NaN*ones([ni*nj,nt]);
for z = 1:50
    id = (kmt(:)==z);
    tb(id,:) = tzall(id,z,:);
end

tb5 = reshape(tb,ni,nj,nt);
time5 = time;

%% Time
%Put times together

temp1 = cat(3,tb1,tb2);
temp2 = cat(3,temp1,tb3);
temp3 = cat(3,temp2,tb4);
temp_btm = cat(3,temp3,tb5);

time = [time1;time2;time3;time4;time5];

%% Select 1950-2007
yr = (time/365)+1888; %days since 1888-01-01 = first month is Jan 1888
runs = find(yr>=1950);
runs(end)
runs(1)
tb = temp_btm(:,:,runs);

%%
save([fpath 'ocean_cobalt_temp_btm_monthly_1908_2007.mat'],'temp_btm','time',...
    'yr','long_name','units','geolat_t','geolon_t');

save([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb','time',...
    'yr','runs','long_name','units','geolat_t','geolon_t');



