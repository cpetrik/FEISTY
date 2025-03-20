% Make GRD file for FEISTY input from CESM-WACCM 1 degree model
% using fixed depth netcdf

clear 
close all

Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/';

%%
load([Cdir 'ukesm_isimip_depth.mat'])

%% Depth, lat, lon
ncdisp([Cdir 'ukesm1-0-ll_r1i1p1f2_picontrol_deptho_onedeg_global_fx.nc'])

ncid = netcdf.open([Cdir 'ukesm1-0-ll_r1i1p1f2_picontrol_deptho_onedeg_global_fx.nc'],'NC_NOWRITE');
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

%% all need to be flipped to match forcings
LAT = fliplr(LAT);
LON = fliplr(LON);
deptho = fliplr(deptho);
lmask = fliplr(lmask);

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,deptho)
title('CMIP')

%%
wID = find(~isnan(deptho(:))); 
nID = length(wID); %41313 (Ndet = 41363) 
%if use daily interp, 50 grid cells without depth info

LID = find(lmask(:)==1);

eq1 = (wID==LID); 
sum(eq1)

%% USE BOTTOM POC TO SET WID AND NID 
% (because that was used in daily interp files)
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/hist/';
load([fpath 'cesm2_hist_det_monthly_1850_2014.mat'],'det');
det(det > 1.0e19) = nan;
test4 = squeeze(double(det(:,:,70)));
figure
pcolor(test4); shading flat
WID = find(~isnan(test4(:)));
NID = length(WID);

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = deptho(ID);
GRD.lmask = lmask(ID);

%% Save needed variables
save([Cdir 'gridspec_cesm2_cmip6_2300.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'Data_grid_cesm2_cmip6_2300.mat'],'GRD');
