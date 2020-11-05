% Compare timeseries of GFDL and IPSL temp and zoop
% Just use surface temp and zmeso vint as indicators

clear all
close all

ipath='/Volumes/MIP/Fish-MIP/CMIP6/IPSL/';
gpath='/Volumes/MIP/Fish-MIP/CMIP6/GFDL/';

% ============================ GFDL ================================
%% Hist
%Meso Zoop vint 
ncid = netcdf.open([gpath 'hist/gfdl-esm4_r1i1p1f1_historical_zmeso-vint_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Time
yr = ((time+1)/12)+1601;
runs = find(yr>1950 & yr<=2015);
test = (length(time)-779):length(time);

% Get subset of time
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs(1)-1], [360 180 length(runs)]);
end
netcdf.close(ncid);
zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;

% time mean
ghzm = squeeze(nanmean(nanmean(zmeso_vint,1),2));

% tos
ncid = netcdf.open([gpath 'hist/gfdl-esm4_r1i1p1f1_historical_tos_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get subset of time
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
%
tos = tos(:,:,runs);
tos(tos >= 1.00e+20) = NaN;

% time mean
ghts = squeeze(nanmean(nanmean(tos,1),2));

clear tos zmeso_vint

%% Pre
%Meso Zoop vint 
ncid = netcdf.open([gpath 'preindust/gfdl-esm4_r1i1p1f1_picontrol_zmeso-vint_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Time
yr = ((time+1)/12)+1600;
runs = find(yr>1950 & yr<=2100);
test = (length(time)-1811):length(time);
t2 = yr(test);

% Get subset of time
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs(1)-1], [360 180 length(runs)]);
end
netcdf.close(ncid);
zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;

% time mean
gpzm = squeeze(nanmean(nanmean(zmeso_vint,1),2));

% tos
ncid = netcdf.open([gpath 'preindust/gfdl-esm4_r1i1p1f1_picontrol_tos_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get subset of time
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
%
tos = tos(:,:,runs);
tos(tos >= 1.00e+20) = NaN;

% time mean
gpts = squeeze(nanmean(nanmean(tos,1),2));

clear tos zmeso_vint

% ============================ IPSL ================================
%% Hist
%Meso Zoop vint 
ncid = netcdf.open([ipath 'hist/ipsl-cm6a-lr_r1i1p1f1_historical_zmeso-vint_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Time
yr = ((time+1)/12)+1601;
runs = find(yr>1950 & yr<=2015);
test = (length(time)-779):length(time);

% Get subset of time
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs(1)-1], [360 180 length(runs)]);
end
netcdf.close(ncid);
zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;

% time mean
ihzm = squeeze(nanmean(nanmean(zmeso_vint,1),2));

% tos
ncdisp([ipath 'hist/ipsl-cm6a-lr_r1i1p1f1_historical_tos_onedeg_global_monthly_1850_2014.nc'])
ncid = netcdf.open([ipath 'hist/ipsl-cm6a-lr_r1i1p1f1_historical_tos_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get subset of time
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
%
tos = tos(:,:,runs);
tos(tos >= 1.00e+20) = NaN;

% time mean
ihts = squeeze(nanmean(nanmean(tos,1),2));

clear tos zmeso_vint

%% Pre
%Meso Zoop vint 
ncid = netcdf.open([ipath 'preindust/ipsl-cm6a-lr_r1i1p1f1_picontrol_zmeso-vint_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Time
yr = ((time+1)/12)+1600;
runs = find(yr>1950 & yr<=2100);
test = (length(time)-1811):length(time);
t2 = yr(test);

% Get subset of time
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs(1)-1], [360 180 length(runs)]);
end
netcdf.close(ncid);
zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;

% time mean
ipzm = squeeze(nanmean(nanmean(zmeso_vint,1),2));

% tos
ncid = netcdf.open([ipath 'preindust/ipsl-cm6a-lr_r1i1p1f1_picontrol_tos_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get subset of time
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
%
tos = tos(:,:,runs);
tos(tos >= 1.00e+20) = NaN;

% time mean
ipts = squeeze(nanmean(nanmean(tos,1),2));

clear tos zmeso_vint

%% SSP 126
%Meso Zoop vint 
ncid = netcdf.open([ipath 'ssp126/ipsl-cm6a-lr_r1i1p1f1_ssp126_zmeso-vint_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

% Time
yr = ((time+1)/12)+1601;
runs = find(yr>1950 & yr<=2015);
test = (length(time)-779):length(time);

% Get subset of time
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    zmeso_vint = netcdf.getVar(ncid,n-1, [0,0,runs(1)-1], [360 180 length(runs)]);
end
netcdf.close(ncid);
zmeso_vint(zmeso_vint >= 1.00e+20) = NaN;

% time mean
ihzm = squeeze(nanmean(nanmean(zmeso_vint,1),2));

% tos
ncid = netcdf.open([ipath 'ssp126/ipsl-cm6a-lr_r1i1p1f1_historical_tos_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get subset of time
for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);
%
tos = tos(:,:,runs);
tos(tos >= 1.00e+20) = NaN;

% time mean
ihts = squeeze(nanmean(nanmean(tos,1),2));

clear tos zmeso_vint

