%Get CMIP6 time for scenarios

clear
close all

%% 126
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp126/';

time_units = 'months since 1601-1-1';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_tob_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Time
t1= time;
yr1 = ((time)/12)+1601;

%%
clear time 

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp126_tob_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Time
t2 = time;
yr2 = ((time)/12)+1601;

%%
clear time 

time = [t1;t2];
yr = [yr1;yr2];

%%
save([fpath 'ipsl_ssp126_time_monthly_2015_2300.mat'],'yr',...
    'time_units','time');

%% 585
clear
close all

fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp585/';

time_units = 'months since 1601-1-1';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_tob_onedeg_global_monthly_2015_2100.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Time
t1= time;
yr1 = ((time)/12)+1601;

%%
clear time 

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp585_tob_60arcmin_global_monthly_2101_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

% Time
t2 = time;
yr2 = ((time)/12)+1601;

%%
clear time 

time = [t1;t2];
yr = [yr1;yr2];

%%
save([fpath 'ipsl_ssp585_time_monthly_2015_2300.mat'],'yr',...
    'time_units','time');

%% 534

clear
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/ssp534over/';

time_units = 'months since 1601-1-1';

%%
ncid = netcdf.open([fpath 'ipsl-cm6a-lr_r1i1p1f1_ssp534-over_thetao-bot_60arcmin_global_monthly_2040_2300.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = ((time)/12)+1601;

%%
save([fpath 'ipsl_ssp534-over_time_monthly_2040_2300.mat'],'yr',...
    'time_units','time');





