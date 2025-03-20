% Make GRD file for FEISTY input from CESM-WACCM 1 degree model
% using fixed depth netcdf

clear 
close all

Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/';

%% Depth, lat, lon
ncdisp([Cdir 'ukesm1-0-ll_deptho_60arcmin_global_fixed.nc'])

ncid = netcdf.open([Cdir 'ukesm1-0-ll_deptho_60arcmin_global_fixed.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%seafloor depths
deptho(deptho >= 1.00e+20) = NaN;
deptho = double(deptho);

lat = double(lat);
lon = double(lon);

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);

%% Land mask
lmask = deptho;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;

%% new 2300 files oriented with 
%  Nhem on left, Shem on right, Russia on top, US on bottom

%% check orientation
figure(4)
subplot(2,2,1)
pcolor(deptho); shading flat
title('dep')

subplot(2,2,3)
pcolor(LAT); shading flat
title('lat')
colorbar

subplot(2,2,4)
pcolor(LON); shading flat
title('lon')
colorbar

subplot(2,2,2)
pcolor(lmask); shading flat
title('lmask')

%% all correct

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)
title('CMIP')

%%
WID = find(~isnan(deptho(:))); 
NID = length(WID); %41363 (Ndet = 41363) 

LID = find(lmask(:)==1);

eq1 = (WID==LID); 
sum(eq1)

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = deptho(ID);
GRD.lmask = lmask(ID);

%% Save needed variables
save([Cdir 'gridspec_ukesm_cmip6_2300.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'Data_grid_ukesm_cmip6_2300.mat'],'GRD');
