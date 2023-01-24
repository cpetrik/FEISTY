% Read GFDL 1 deg netcdfs
% obsclim
% time

clear 
close all

%fpath='/Volumes/MIP/Fish-MIP/Phase3/OneDeg/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

%% one file
ncdisp([fpath 'gfdl-mom6-cobalt2_obsclim_tob_onedeg_global_monthly_1961_2010.nc'])

%%
% time
% Size:       600x1
% Dimensions: time
% Datatype:   double
% Attributes:
time_standard_name = 'time';
time_long_name     = 'time';
time_units         = 'months since 1901-1-1 00:00:00';
calendar           = '360_day';
time_axis          = 'T';

%%
ncid = netcdf.open([fpath 'gfdl-mom6-cobalt2_obsclim_tob_onedeg_global_monthly_1961_2010.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
%Change to time since 1841

% Hist
yr_hist = 1901 + ((time+1)/12);
time_hist = (yr_hist - 1841) * 12 - 1;

% PI/trans time 
yr_trans = 1841+(1/12):(1/12):1961;
time_trans = (yr_trans - 1841) * 12 - 1;

% spinup last 10 yrs
yr_spin = 1831+(1/12):(1/12):1841;
time_spin = (yr_spin - 1841) * 12 - 1;

%%
save([fpath 'gfdl-mom6-cobalt2_obsclim_time_onedeg_global_monthly_1831_2010.mat'],...
    'time_long_name','time_standard_name','time_units','calendar','time_axis',...
    'lat','lon','time_hist','yr_hist','time_trans','yr_trans','time_spin','yr_spin');
