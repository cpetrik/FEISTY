% Make GRD file for FEISTY input from IPSL 1 degree model
% using fixed depth netcdf

clear 
%close all

Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/hist/';

%%
load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/gridspec_ipsl_cmip6.mat','deptho')

%% Depth, lat, lon
ncdisp([Cdir 'IPSL/ipsl-cm6a-lr_r1i1p1f1_picontrol_deptho_onedeg_global_fixed.nc'])

ncid = netcdf.open([Cdir 'IPSL/ipsl-cm6a-lr_r1i1p1f1_picontrol_deptho_onedeg_global_fixed.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%seafloor depths
deptho(deptho >= 1.00e+20) = NaN;
deptho = double(deptho);

%Grid of lat & lon
[LAT,LON] = meshgrid(lat, lon);

%% check orientation
figure
pcolor(deptho)

figure
pcolor(LAT)

LAT = fliplr(LAT);
deptho = fliplr(deptho);

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
NID = length(WID); %41383

% Land mask
lmask = deptho;
lmask(~isnan(lmask)) = 1;
lmask(isnan(lmask)) = 0;
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
save([Cdir 'IPSL/gridspec_ipsl_cmip6.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'IPSL/Data_grid_ipsl.mat'],'GRD');
