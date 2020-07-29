% Read Fish-MIP netcdfs
% GFDL preindust 1601-2100

clear all
close all

fpath='/Volumes/FEISTY/Fish-MIP/CMIP6/GFDL/preindust/';

%% Bottom temp
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_picontrol_thetao-bot_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%thetao: 360 x 180 x 6000 (500 yrs)
%time: 6000
%NaNs = 1.0000e+20

%% Bottom det
ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_picontrol_expc-bot_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%expc: 360 x 180 x 6000 (500 yrs)
%time: 6000
%NaNs = 1.0000e+20

%% SST
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_picontrol_tos_onedeg_global_monthly_1601_2100.nc'])

ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_picontrol_tos_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%tos: 360 x 180 x 6000 (500 yrs)
%time: 6000
%NaNs = 1.0000e+20

%% MesoZ zint
% ncdisp([fpath 'gfdl-esm4_r1i1p1f1_picontrol_zmeso-vint_onedeg_global_monthly_1601_2100.nc'])
%'Mole Content of Mesozooplankton expressed as Carbon in sea water'
%units    = 'mol m-2'
%comment  = 'vertically integrated over all ocean levels by ISIMIP data management team'

ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_picontrol_zmeso-vint_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    if i==nvars
        zmeso_vint=netcdf.getVar(ncid,i-1);
        zmeso_vint(zmeso_vint == 1.0000e+20) = NaN;
    else
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end
end
netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%zmeso-vint: 360 x 180 x 6000 (500 yrs) MATLAB won't read in with dash "-"  
%time: 6000
%NaNs = 1.0000e+20

%% MicroZ zint
ncdisp([fpath 'gfdl-esm4_r1i1p1f1_picontrol_zmicro-vint_onedeg_global_monthly_1601_2100.nc'])
%'Mole Content of Microzooplankton expressed as Carbon in sea water'
%units    = 'mol m-2'
%comment  = 'vertically integrated over all ocean levels by ISIMIP data management team'

ncid = netcdf.open([fpath 'gfdl-esm4_r1i1p1f1_picontrol_zmicro-vint_onedeg_global_monthly_1601_2100.nc'],'NC_NOWRITE');
%fname='mom6_preindust_gridspec';

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    if i==nvars
        zmicro_vint=netcdf.getVar(ncid,i-1);
        zmicro_vint(zmicro_vint == 1.000000020040877e+20) = NaN;
    else
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
    end
end
netcdf.close(ncid);

% Vars
%lat: 180
%lon: 360
%zmicro-vint: 360 x 180 x 6000 (500 yrs) MATLAB won't read in with dash "-"  
%time: 6000
%NaNs = 1.0000e+20

     



