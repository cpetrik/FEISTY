% Calc thkcello for MOM6-NEP

clear 
close all

%% 
fpath='/project/Feisty/GCM_Data/MOM6-NEP10/';

ncdisp([fpath 'ocean_static.nc'])

ncdisp([fpath 'thetao.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc'])

%%
%ocean_static
% areacello     var2
% Size:       775x845
% Dimensions: ih,jh
% Datatype:   single
% Attributes:
% _FillValue    = 1.000000020040877e+20
% missing_value = 1.000000020040877e+20
% units         = 'm2'
% long_name     = 'Ocean Grid-Cell Area'
% cell_methods  = 'area:sum jh:sum ih:sum time: point'
% standard_name = 'cell_area'

%thetao
% volcello      var8
% Size:       775x845x52x372
% Dimensions: xh,yh,z_l,time
% Datatype:   single
% Attributes:
% _FillValue    = 1.000000020040877e+20
% missing_value = 1.000000020040877e+20
% units         = 'm3'
% long_name     = 'Ocean grid-cell volume'
% cell_methods  = 'area:sum z_l:sum yh:sum xh:sum time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% standard_name = 'ocean_volume'

% time
time_units    = 'days since 1993-01-01 00:00:00';
calendar_type = 'GREGORIAN';
calendar      = 'gregorian';

%% Area
ncid = netcdf.open([fpath 'ocean_static.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 2
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

areacello(areacello>1e20) = nan;

%% time & z_l
ncid = netcdf.open([fpath 'thetao.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% lat & lon
for i = 4:5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
%time
for i = 8
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end

%find cell thicknesses to get top 100m
for i = 12
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end

z_l(z_l>1e20) = nan;

%% just top 100 of volcello
z100 = find(z_l <= 100);

ni = length(ih);
nj = length(jh);
nz = length(z_l);
nt = length(time);

for i = 10
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj length(z100) nt]);']);
    eval([ varname '(' varname ' == 1.0000000e+20) = NaN;']);
end
netcdf.close(ncid);

volcello(volcello>1e20) = nan;

%% Thk = volcello./areacello

areamat = repmat(areacello,1,1,length(z100),nt);

%%
thkcello = volcello ./ areacello;

thk_units = 'm';
zl = z_l(z100);

save([fpath 'thkcello.nep.full.hcast.monthly.raw.r20250912.199301-202506.mat'],...
    'thkcello','thk_units','ih','jh','zl','time', '-v7.3');
