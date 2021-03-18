% Read GFDL ESM4 Preindust netcdfs
% 3D temp

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/ESM2M_PI/Temp3D/';

%%
ncdisp([fpath 'ocean.100101-102012.temp.nc'])

% Dimensions:
% xt_ocean       = 360
% yt_ocean       = 200
% time           = 240   (UNLIMITED)
% nv             = 2
% xu_ocean       = 360
% yu_ocean       = 200
% st_ocean       = 50
% st_edges_ocean = 51
% time units     = 'days since 0001-01-01 00:00:00'

% temp
% Size:       360x200x50x240
% Dimensions: xt_ocean,yt_ocean,st_ocean,time
% Datatype:   single
% Attributes:
long_name     = 'Potential temperature';
units         = 'degrees K';
% valid_range   = [-10  500]
% missing_value = -1.000000020040877e+20
% _FillValue    = -1.000000020040877e+20
% cell_methods  = 'time: mean'
% time_avg_info = 'average_T1,average_T2,average_DT'
% coordinates   = 'geolon_t geolat_t'
% standard_name = 'sea_water_potential_temperature'

%% Dates
s20 = 1001:20:1381;
e20 = 1020:20:1400;

%% Cycle through each 20 yrs
temp_100 = [];
time_all = [];
for q = 1:length(s20)
    
    tspan = [num2str(s20(q)),'01-',num2str(e20(q)),'12'];
    
    %%
    ncid = netcdf.open([fpath 'ocean.',tspan,'.temp.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == -1.0e+10) = NaN;']);
    end
    netcdf.close(ncid);
    
    %% Top 100 m
    z100 = find(st_ocean <= 100);
    temp(temp <= -1e19) = NaN;
    tp = temp(:,:,z100,:);
    
    %% Mean top 100 m
    %constant thickness of 10 m in top 200 m
    tp_100 = squeeze(nanmean(tp,3));
    
    %% save
    save([fpath 'ocean.',tspan,'.temp_100.mat'],'tp_100','time',...
        'geolon_t','geolat_t','long_name','units','st_ocean','z100')
    
    clear temp tp tp_100 time
    
end

