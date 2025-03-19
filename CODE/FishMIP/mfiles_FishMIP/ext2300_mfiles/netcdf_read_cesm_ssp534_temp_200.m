% Read CMIP6 netcdfs
% CESM SSP534-over
% Mean top 0-150 m

clear 
close all

%fpath='/project/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/ssp534over/';
tpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';

%% zoo zall
ncdisp([fpath 'cesm2-waccm_r1i1p1f1_ssp534-over_thetao_60arcmin_global_monthly_2040_2299.nc'])

%%
standard_name = 'sea_water_potential_temperature';
long_name     = 'Sea Water Potential Temperature';
units         = 'degC';
missing_value = 1.000000020040877e+20;
%Size:       320x384x15x1980
%Dimensions: i,j,lev,time
%time units = 'months since 1601-1-1 00:00:00'

%% 
ncid = netcdf.open([fpath 'cesm2-waccm_r1i1p1f1_ssp534-over_thetao_60arcmin_global_monthly_2040_2299.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get all other vars 1st
for n = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

%% Get subsets of thetao so smaller memory?
yr = ((time+1)/12)+1601;
%z200 = find(lev <= 150e2); %lev in cm
z200 = 1:15; % only top 15 levels saved in plankton

for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    thetao = netcdf.getVar(ncid,n-1, [0,0,z200(1)-1,0], [360 180 length(z200) length(time)]);
end
netcdf.close(ncid);
thetao(thetao >= 1.00e+20) = NaN;

%% Subset of thkcello #lev_bnds in m (not cm like lev)
[ni,nj,~,nt] = size(thetao);
thkcello = (lev_bnds(2,z200) - lev_bnds(1,z200))';
thkcello = repmat(thkcello,1,ni,nj,nt);
thkcello = permute(thkcello,[2 3 1 4]);

%% Mean top 150 m
temp_150 = squeeze(sum((thetao.*thkcello),3,'omitnan')) ./ squeeze(sum((thkcello),3,'omitnan'));

%% Clear 
clear thetao thkcello

save([fpath 'cesm2_ssp534-over_temp_150_monthly_2040_2299.mat'],'temp_150',...
    'yr','long_name','standard_name','units','lat','lon',...
    'lev','lev_bnds','time','z200');
