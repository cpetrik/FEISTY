% Make GRD file for FEISTY input from GFDL 1/4 degree 
% Reanalysis hindcast model

clear all
close all

Cdir = '/Volumes/MIP/Fish-MIP/Phase3/';

%% ncinfo
ncdisp([Cdir 'gfdl-mom6-cobalt2_obsclim_deptho_15arcmin_global_fixed.nc'])

%%
% Global Attributes:
%            CDI                      = 'Climate Data Interface version 1.9.8 (https://mpimet.mpg.de/cdi)'
%            Conventions              = 'CF-1.6'
%            filename                 = '19880101.ocean_static.nc'
%            title                    = 'OM4_SIS2_COBALT_initFromSrc_grid'
%            grid_type                = 'regular'
%            grid_tile                = 'N/A'
%            cdo_openmp_thread_number = 4
%            isimip_id                = '40c5153d-8512-4afc-8c96-4149e92b366d'
% Dimensions:
%            lon = 1440
%            lat = 720
% Variables:
%     deptho
%            Size:       1440x720
%            Dimensions: lon,lat
standard_name = 'sea_floor_depth_below_geoid';
long_name     = 'Sea Floor Depth Below Geoid';
units         = 'm';
FillValue    = -9.999999790214768e+33;
missing_value = -9.999999790214768e+33;
GFDL_variable = 'deptho';

%% Depth, Lat & Lon
ncid = netcdf.open([Cdir 'gfdl-mom6-cobalt2_obsclim_deptho_15arcmin_global_fixed.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%%
deptho(deptho <= -1.00e+33) = NaN;
deptho = double(deptho);

%% Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);

%%
figure
pcolor(deptho); shading flat
colorbar

figure
pcolor(LAT); shading flat
colorbar

figure
pcolor(LON); shading flat
colorbar

%%
close all

%% Land mask
lmask = deptho;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

LID = find(lmask(:)==1);
WID = find(~isnan(deptho));  % spatial index of water cells
NID = length(WID);

eq1 = (WID==LID);
sum(eq1)

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = deptho(ID);
GRD.lmask = lmask(ID);

%% Save needed variables
save([Cdir 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin_.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin_.mat'],'GRD');
