% Make GRD file for FEISTY input from NEMURO CCE model
% using IPSL downscaled projections

clear 
close all

%% IPSL
Cdir = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
load([Cdir 'feisty_ipsl_gridspec.mat'],'LAT')
load([Cdir 'Data_grid_nemuro_ipsl.mat'],'GRD')

% Area calc
%It is a regular 0.1x0.1 deg grid, so you can easily calculate the area (in km2) of each grid cell as:
%area = dx * dy = (0.1*111.1*cos(lat)) * (0.1*111.1) = 123.43*cos(lat)
%convert latitude to radian first? 
%Latitude is between 0-90deg (0-pi/2), so the cos should be positive.
rLAT = LAT.*pi./180;
Area = abs(123.43*cos(rLAT));

% Retain only water cells
ID = GRD.ID;
GRD.AREA  = Area(ID);

%% Save needed variables
save([Cdir 'feisty_ipsl_gridspec.mat'],'Area','-append');
save([Cdir 'Data_grid_nemuro_ipsl.mat'],'GRD');

clear rLAT LAT ID Area GRD

%% GFDL
Gdir = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/GFDLdown/';
load([Gdir 'feisty_gfdl_gridspec.mat'],'LAT')
load([Gdir 'Data_grid_nemuro_gfdl.mat'],'GRD')

% Area calc
rLAT = LAT.*pi./180;
Area = abs(123.43*cos(rLAT));

% Retain only water cells
ID = GRD.ID;
GRD.AREA  = Area(ID);

% Save needed variables
save([Gdir 'feisty_gfdl_gridspec.mat'],'Area','-append');
save([Gdir 'Data_grid_nemuro_gfdl.mat'],'GRD');

clear rLAT LAT ID Area GRD

%% HAD
Hdir = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/HADdown/';
load([Hdir 'feisty_hadley_gridspec.mat'],'LAT')
load([Hdir 'Data_grid_nemuro_hadley.mat'],'GRD')

% Area calc
rLAT = LAT.*pi./180;
Area = abs(123.43*cos(rLAT));

% Retain only water cells
ID = GRD.ID;
GRD.AREA  = Area(ID);

% Save needed variables
save([Hdir 'feisty_hadley_gridspec.mat'],'Area','-append');
save([Hdir 'Data_grid_nemuro_hadley.mat'],'GRD');

clear rLAT LAT ID Area GRD
