% Read CMIP6 netcdfs
% CESM deptho

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp126/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp126/';

%% deptho
ncdisp([fpath 'cesm2-waccm_r1i1p1f1_ssp126_deptho_60arcmin_global_fx.nc'])

%%
units         = 'm'; %? or cm?
missing_value = 1.000000020040877e+20;
%Size:       320x384x15x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1 00:00:00'

%% 
ncid = netcdf.open([fpath 'cesm2-waccm_r1i1p1f1_ssp126_deptho_60arcmin_global_fx.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

deptho(deptho >= 1.00e+20) = NaN;

save([fpath 'cesm2-waccm_r1i1p1f1_ssp126_deptho_60arcmin_global_fx.mat']);
