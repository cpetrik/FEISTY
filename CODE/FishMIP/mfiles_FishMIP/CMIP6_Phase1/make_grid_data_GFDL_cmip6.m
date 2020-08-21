% Make GRD file for FEISTY input from CESM 1 degree model

clear all
close all

Cdir = '/Volumes/FEISTY/Fish-MIP/CMIP6/';

%% Lat & Lon
ncid = netcdf.open([Cdir 'GFDL/hist/gfdl-esm4_r1i1p1f1_historical_thetao-bot_onedeg_global_monthly_1850_2014.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);
thetao(thetao >= 1.00e+20) = NaN;

%Land mask
mask = squeeze(thetao(:,:,1));
lmask = double(mask);
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);

LID = find(lmask(:)==1);
WID = find(~isnan(thetao(:,:,1)));  % spatial index of water cells
NID = length(WID);

eq1 = (WID==LID);
sum(eq1)

%% Depth, lat, lon
ncid = netcdf.open([Cdir 'GFDL/gfdl-esm4_r1i1p1f1_picontrol_deptho_onedeg_global_fx.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%%
deptho(deptho >= 1.00e+20) = NaN;
deptho = double(deptho);

%%
DID = find(~isnan(deptho(:))); 
IDN = length(DID);

eq2 = (WID==DID);
sum(eq2)

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
save([Cdir 'GFDL/gridspec_gfdl_cmip6.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'GFDL/Data_grid_gfdl.mat'],'GRD');
