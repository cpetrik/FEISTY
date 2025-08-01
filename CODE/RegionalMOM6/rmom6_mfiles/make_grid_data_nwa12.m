% Make grid data for NWA12
% TP and vel have fewer ocean cells than BGC
% Comp to wet?

clear
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NWA12/';

%%
load([fpath 'nwa_raw_ocean_static_gridspec.mat'])

%% index of water cells
WID = find((wet(:)==1));
WIDc = find((wet_c(:)==1));
WIDu = find((wet_u(:)==1));
WIDv = find((wet_v(:)==1));

WIDtu = intersect(WID,WIDu);
WIDtuv = intersect(WIDv,WIDtu);

%%
NID = length(WID);

GRD.ID = WID;
GRD.NID = NID;
GRD.Lat = double(geolat(WID));
GRD.Lon = double(geolon(WID));
GRD.depth = double(deptho(WID));
GRD.area = double(areacello(WID));

save([fpath 'Data_grid_mom6_nwa12.mat'], 'GRD');

%%
[ni,nj] = size(geolon);
mask = zeros(ni,nj);
mask(WID) = ones;

GRD2.area = double(areacello);
GRD2.dx = double(dxt);
GRD2.dy = double(dyt);
GRD2.lat = double(geolat);
GRD2.lon = double(geolon);
GRD2.mask = double(mask);

save([fpath 'Data_grid2D_mom6_nwa12.mat'], 'GRD2');

%%
% load([fpath 'Data_BGC_grid_mom6_nwa12.mat'], 'GRD');
% GRDv0 = GRD;
% clear GRD
% 
% %% Some grid cells with wet=1 have nan velocities
% % UID = find(~isnan(ut_100(:,:,1)));
% % VID = find(~isnan(vt_100(:,:,1)));
% % UVID = intersect(UID,VID);
% 
% GRD.UVID = GRDv0.UVID;
% UVID = GRDv0.UVID;
% NID = length(UVID);
% 


