% Read MOM6 ocean static netcdfs

clear 
close all

%% NEP
ppath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';

% 
ncdisp([ppath 'ocean_static.nc'])

%% NWA
apath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';

% 
ncdisp([apath 'ocean_static.nc'])