% Read GFDL ESM2M Preindust netcdfs
% all vars for FEISTY

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/ESM2M_PI/';
tpath='/Volumes/MIP/GCM_DATA/ESM2M_PI/Temp3D/';

%%
ncdisp([fpath 'ocean_cobalt_btm.100101-102012.fndet_btm.nc'])

% Dimensions:
% xu_ocean = 360
% yu_ocean = 200
% time     = 240   (UNLIMITED)
% nv       = 2
% xt_ocean = 360
% yt_ocean = 200
% time:units = 'days since 0001-01-01 00:00:00'

% fndet_btm
% Size:       360x200x240
% Dimensions: xt_ocean,yt_ocean,time
% Datatype:   single
% Attributes:
det_long_name     = 'ndet sinking flux to bottom';
det_units         = 'mol m-2 s-1';
% missing_value = -10000000000
% _FillValue    = -10000000000
% cell_methods  = 'time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% coordinates   = 'geolon_t geolat_t'

%%
% btm_temp
tb_long_name     = 'Bottom Temperature';
tb_units         = 'deg C';

% nmdz_100
mz_long_name     = 'Medium zooplankton nitrogen biomass in upper 100m';
mz_units         = 'mol m-2';

% nlgz_100
lz_long_name     = 'Large zooplankton nitrogen biomass in upper 100m';
lz_units         = 'mol m-2';

% jhploss_nmdz_100
mzhp_long_name     = 'Medium zooplankton nitrogen loss to higher preds. integral in upper 100m';
mzhp_units         = 'mol m-2 s-1';

% jhploss_nlgz_100
lzhp_long_name     = 'Large zooplankton nitrogen loss to higher preds. integral in upper 100m';
lzhp_units         = 'mol m-2 s-1';

%% Dates
%s20 = 1001:20:1481;
%e20 = 1020:20:1500; only 400 yrs of 3d temp

s20 = 1001:20:1381;
e20 = 1020:20:1400;

%% Cycle through each 20 yrs
time_all = [];
Tp = [];
Tb = [];
mz = [];
lz = [];
mzhp = [];
lzhp = [];
det = [];
for q = 1:length(s20)
    
    tspan = [num2str(s20(q)),'01-',num2str(e20(q)),'12'];
    
    %%
    load([tpath 'ocean.',tspan,'.temp_100.mat']);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_btm.',tspan,'.fndet_btm.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_btm.',tspan,'.btm_temp.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_biomass_100.',tspan,'.nmdz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_biomass_100.',tspan,'.nlgz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_miscflux_100.',tspan,'.jhploss_nmdz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    ncid = netcdf.open([fpath 'ocean_cobalt_miscflux_100.',tspan,'.jhploss_nlgz_100.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    %% NaNs on land cells
    fndet_btm(fndet_btm <= -1e9) = NaN;
    btm_temp(btm_temp <= -1e9) = NaN;
    nmdz_100(nmdz_100 <= -1e9) = NaN;
    nlgz_100(nlgz_100 <= -1e9) = NaN;
    jhploss_nmdz_100(jhploss_nmdz_100 <= -1e9) = NaN;
    jhploss_nlgz_100(jhploss_nlgz_100 <= -1e9) = NaN;
    
    %% Convert units after interpolating
    %det btm: mol N m-2 s-1
    %zoo: integrated to mol N m-2
    %zoo loss: mol N m-2 s-1
    %tp: K
    %tb: degC
    
    tp_100 = double(tp_100) - 273;
    btm_temp = double(btm_temp);
    nmdz_100 = double(nmdz_100) * (106.0/16.0) * 12.01 * 9.0;
    nlgz_100 = double(nlgz_100) * (106.0/16.0) * 12.01 * 9.0;
    jhploss_nmdz_100 = double(jhploss_nmdz_100) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
    jhploss_nlgz_100 = double(jhploss_nlgz_100) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
    fndet_btm = double(fndet_btm) * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
    
    %% Concat
    time_all = [time_all; time];
    Tp = cat(3,Tp,tp_100);
    Tb = cat(3,Tb,btm_temp);
    mz = cat(3,mz,nmdz_100);
    lz = cat(3,lz,nlgz_100);
    mzhp = cat(3,mzhp,jhploss_nmdz_100);
    lzhp = cat(3,lzhp,jhploss_nlgz_100);
    det = cat(3,det,fndet_btm);
    
end

%% Negative biomass or mortality loss from interp
mz(mz<0) = 0.0;
lz(lz<0) = 0.0;
mzhp(mzhp<0) = 0.0;
lzhp(lzhp<0) = 0.0;
det(det<0) = 0.0;

%% Save
% save([fpath 'ocean_cobalt_1001_1400.mat'],'time_all','Tp','Tb','mz','lz',...
%     'mzhp','lzhp','det');

save([fpath 'ocean_cobalt_temp_1001_1400.mat'],'time_all','Tp','Tb','-v7.3');

save([fpath 'ocean_cobalt_det_1001_1400.mat'],'time_all','det','-v7.3');

save([fpath 'ocean_cobalt_zoo_1001_1400.mat'],'time_all','mz','lz','-v7.3');

save([fpath 'ocean_cobalt_hploss_1001_1400.mat'],'time_all','mzhp','lzhp','-v7.3');







