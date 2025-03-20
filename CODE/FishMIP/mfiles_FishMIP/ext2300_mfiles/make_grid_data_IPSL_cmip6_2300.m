% Make GRD file for FEISTY input from IPSL 1 degree model
% using fixed depth netcdf

clear 
%close all

Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/IPSL/';

%%
load([Cdir 'gridspec_ipsl_cmip6.mat'])

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
WID = find(~isnan(deptho(:))); 
NID = length(WID); %41383

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
save([Cdir 'gridspec_ipsl_cmip6_2300.mat'],'deptho','LAT','LON','lmask');
save([Cdir 'Data_grid_ipsl_cmip6_2300.mat'],'GRD');
